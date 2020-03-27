#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <queue>
#include <atomic>
#include <vector>

#include <cuda_runtime.h>

#include "converse.h"
#include "hapi.h"
#include "hapi_impl.h"
#include "gpumanager.h"
#ifdef HAPI_NVTX_PROFILE
#include "hapi_nvtx.h"
#endif

#if defined HAPI_TRACE || defined HAPI_INSTRUMENT_WRS
extern "C" double CmiWallTimer();
#endif

#ifdef HAPI_TRACE
#define QUEUE_SIZE_INIT 128
extern "C" int traceRegisterUserEvent(const char* x, int e);
extern "C" void traceUserBracketEvent(int e, double beginT, double endT);

typedef struct gpuEventTimer {
  int stage;
  double cmi_start_time;
  double cmi_end_time;
  int event_type;
  const char* trace_name;
} gpuEventTimer;
#endif

static void createPool(int *nbuffers, int n_slots, std::vector<BufferPool> &pools);
static void releasePool(std::vector<BufferPool> &pools);

// Event stages used for profiling.
enum WorkRequestStage{
  DataSetup        = 1,
  KernelExecution  = 2,
  DataCleanup      = 3
};

enum ProfilingStage{
  GpuMemSetup   = 8800,
  GpuKernelExec = 8801,
  GpuMemCleanup = 8802
};

#ifndef HAPI_CUDA_CALLBACK
typedef struct hapiEvent {
  cudaEvent_t event;
  void* cb;
  void* cb_msg;
  hapiWorkRequest* wr; // if this is not NULL, buffers and request itself are deallocated

  hapiEvent(cudaEvent_t event_, void* cb_, void* cb_msg_, hapiWorkRequest* wr_ = NULL)
            : event(event_), cb(cb_), cb_msg(cb_msg_), wr(wr_) {}
} hapiEvent;

CpvDeclare(std::queue<hapiEvent>, hapi_event_queue);
#endif
CpvDeclare(int, n_hapi_events);

// Used to invoke user's Charm++ callback function
void (*hapiInvokeCallback)(void*, void*) = NULL;

// Functions used to support quiescence detection.
void (*hapiQdCreate)(int) = NULL;
void (*hapiQdProcess)(int) = NULL;

#define MAX_PINNED_REQ 64
#define MAX_DELAYED_FREE_REQS 64

// Declare GPU Manager as a process-shared object.
CsvDeclare(GPUManager, gpu_manager);

CpvDeclare(int, my_device); // GPU device that this thread is mapped to
CpvDeclare(bool, device_rep); // Is this PE a device representative thread? (1 per device)

// Initialize per-process variables
void hapiInitCsv() {
  // Create and initialize GPU Manager object
  CsvInitialize(GPUManager, gpu_manager);
  CsvAccess(gpu_manager).init();
}

// Initialize per-PE variables
void hapiInitCpv() {
  // HAPI event-related
#ifndef HAPI_CUDA_CALLBACK
  CpvInitialize(std::queue<hapiEvent>, hapi_event_queue);
#endif
  CpvInitialize(int, n_hapi_events);
  CpvAccess(n_hapi_events) = 0;

  // Device mapping
  CpvInitialize(int, my_device);
  CpvAccess(my_device) = 0;
  CpvInitialize(bool, device_rep);
  CpvAccess(device_rep) = false;
}

// Clean up per-process data
void hapiExitCsv() {
  // Destroy GPU Manager object
  CsvAccess(gpu_manager).destroy();

  // Release memory pool
  if (CsvAccess(gpu_manager).mempool_initialized_) {
    releasePool(CsvAccess(gpu_manager).mempool_free_bufs_);
  }
}

// Set up PE to GPU mapping, invoked from all PEs
// TODO: Support custom mappings
void hapiMapping(char** argv) {
  Mapping map_type = Mapping::Block; // Default is block mapping
  bool all_gpus = false; // If true, all GPUs are visible to all processes.
                         // Otherwise, only a subset are visible (e.g. with jsrun)
  char* gpumap = NULL;

  // Process +gpumap
  if (CmiGetArgStringDesc(argv, "+gpumap", &gpumap,
        "define pe to gpu device mapping")) {
    if (CmiMyPe() == 0) {
      CmiPrintf("HAPI> PE-GPU mapping: %s\n", gpumap);
    }

    if (strcmp(gpumap, "none") == 0) {
      map_type = Mapping::None;
    }
    else if (strcmp(gpumap, "block") == 0) {
      map_type = Mapping::Block;
    }
    else if (strcmp(gpumap, "roundrobin") == 0) {
      map_type = Mapping::RoundRobin;
    }
    else {
      CmiAbort("Unsupported mapping type: %s, use one of \"none\", \"block\", "
          "\"roundrobin\"", gpumap);
    }
  }

  // Process +allgpus
  if (CmiGetArgFlagDesc(argv, "+allgpus",
        "all GPUs are visible to all processes")) {
    all_gpus = true;
    if (CmiMyPe() == 0) {
      CmiPrintf("HAPI> All GPUs are visible to all processes\n");
    }
  }

  // No mapping specified, user assumes responsibility
  if (map_type == Mapping::None) {
    if (CmiMyPe() == 0) {
      CmiPrintf("HAPI> User should explicitly select devices for PEs/chares\n");
    }
    return;
  }

  CmiAssert(map_type != Mapping::None);

  if (CmiMyRank() == 0) {
    // Count number of GPU devices used by each process
    int visible_device_count;
    hapiCheck(cudaGetDeviceCount(&visible_device_count));
    if (visible_device_count <= 0) {
      CmiAbort("Unable to perform PE-GPU mapping, no GPUs found!");
    }

    int& device_count = CsvAccess(gpu_manager).device_count;
    if (all_gpus) {
      device_count = visible_device_count / (CmiNumNodes() / CmiNumPhysicalNodes());
    } else {
      device_count = visible_device_count;
    }

    // Handle the case where the number of GPUs per process are larger than
    // the number of PEs per process. This is needed because we currently don't
    // support each PE using more than one device.
    if (device_count > CmiNodeSize(CmiMyNode())) {
      if (CmiMyPe() == 0) {
        CmiPrintf("HAPI> Found more GPU devices (%d) than PEs (%d) per process, "
            "limiting to %d device(s) per process\n", device_count,
            CmiNodeSize(CmiMyNode()), CmiNodeSize(CmiMyNode()));
      }
      device_count = CmiNodeSize(CmiMyNode());
    }

    // Count number of PEs per device
    CsvAccess(gpu_manager).pes_per_device = CmiNodeSize(CmiMyNode()) / device_count;

    // Count number of devices on a physical node
    CsvAccess(gpu_manager).device_count_on_physical_node =
      device_count * (CmiNumNodes() / CmiNumPhysicalNodes());
  }

  if (CmiMyPe() == 0) {
    CmiPrintf("HAPI> Config: %d device(s) per process, %d PE(s) per device, %d device(s) per host\n",
        CsvAccess(gpu_manager).device_count, CsvAccess(gpu_manager).pes_per_device,
        CsvAccess(gpu_manager).device_count_on_physical_node);
  }

  CmiNodeBarrier();

  // Perform mapping and set device representative PE
  int my_rank = all_gpus ? CmiPhysicalRank(CmiMyPe()) : CmiMyRank();

  switch (map_type) {
    case Mapping::Block:
      CpvAccess(my_device) = my_rank / CsvAccess(gpu_manager).pes_per_device;
      if (my_rank % CsvAccess(gpu_manager).pes_per_device == 0) CpvAccess(device_rep) = true;
      break;
    case Mapping::RoundRobin:
      CpvAccess(my_device) = my_rank % CsvAccess(gpu_manager).device_count;
      if (my_rank < CsvAccess(gpu_manager).device_count) CpvAccess(device_rep) = true;
      break;
    default:
      CmiAbort("Unsupported mapping type!");
  }

  // Set device for each PE
  hapiCheck(cudaSetDevice(CpvAccess(my_device)));
}

#ifndef HAPI_CUDA_CALLBACK
void recordEvent(cudaStream_t stream, void* cb, void* cb_msg, hapiWorkRequest* wr = NULL) {
  // create CUDA event and insert into stream
  cudaEvent_t ev;
  cudaEventCreateWithFlags(&ev, cudaEventDisableTiming);
  cudaEventRecord(ev, stream);

  hapiEvent hev(ev, cb, cb_msg, wr);

  // push event information in queue
  CpvAccess(hapi_event_queue).push(hev);

  // increase count so that scheduler can poll the queue
  CpvAccess(n_hapi_events)++;
}
#endif

inline static void hapiWorkRequestCleanup(hapiWorkRequest* wr) {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).progress_lock_);
#endif

  // free device buffers
  CsvAccess(gpu_manager).freeBuffers(wr);

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).progress_lock_);
#endif

  // free hapiWorkRequest
  delete wr;
}

#ifdef HAPI_CUDA_CALLBACK
// Invokes user's host-to-device callback.
static void* hostToDeviceCallback(void* arg) {
#ifdef HAPI_NVTX_PROFILE
  NVTXTracer nvtx_range("hostToDeviceCallback", NVTXColor::Asbestos);
#endif
  hapiWorkRequest* wr = *((hapiWorkRequest**)((char*)arg + CmiMsgHeaderSizeBytes + sizeof(int)));
  CmiAssert(hapiInvokeCallback);
  hapiInvokeCallback(wr->host_to_device_cb);

  // inform QD that the host-to-device transfer is complete
  CmiAssert(hapiQdProcess);
  hapiQdProcess(1);

  return NULL;
}

// Invokes user's kernel execution callback.
static void* kernelCallback(void* arg) {
#ifdef HAPI_NVTX_PROFILE
  NVTXTracer nvtx_range("kernelCallback", NVTXColor::Asbestos);
#endif
  hapiWorkRequest* wr = *((hapiWorkRequest**)((char*)arg + CmiMsgHeaderSizeBytes + sizeof(int)));
  CmiAssert(hapiInvokeCallback);
  hapiInvokeCallback(wr->kernel_cb);

  // inform QD that the kernel is complete
  CmiAssert(hapiQdProcess);
  hapiQdProcess(1);

  return NULL;
}

// Frees device buffers and invokes user's device-to-host callback.
// Invoked regardless of the availability of the user's callback.
static void* deviceToHostCallback(void* arg) {
#ifdef HAPI_NVTX_PROFILE
  NVTXTracer nvtx_range("deviceToHostCallback", NVTXColor::Asbestos);
#endif
  hapiWorkRequest* wr = *((hapiWorkRequest**)((char*)arg + CmiMsgHeaderSizeBytes + sizeof(int)));

  // invoke user callback
  if (wr->device_to_host_cb) {
    CmiAssert(hapiInvokeCallback);
    hapiInvokeCallback(wr->device_to_host_cb);
  }

  hapiWorkRequestCleanup(wr);

  // inform QD that device-to-host transfer is complete
  CmiAssert(hapiQdProcess);
  hapiQdProcess(1);

  return NULL;
}

// Used by lightweight HAPI.
static void* lightCallback(void *arg) {
#ifdef HAPI_NVTX_PROFILE
  NVTXTracer nvtx_range("lightCallback", NVTXColor::Asbestos);
#endif

  char* conv_msg_tmp = (char*)arg + CmiMsgHeaderSizeBytes + sizeof(int);
  void* cb = *((void**)conv_msg_tmp);

  // invoke user callback
  if (cb) {
    CmiAssert(hapiInvokeCallback);
    hapiInvokeCallback(cb);
  }

  // notify process to QD
  CmiAssert(hapiQdProcess);
  hapiQdProcess(1);

  return NULL;
}
#endif // HAPI_CUDA_CALLBACK

// Register callback functions. All PEs need to call this.
void hapiRegisterCallbacks() {
#ifdef HAPI_CUDA_CALLBACK
  // FIXME: Potential race condition on assignments, but CmiAssignOnce
  // causes a hang at startup.
  CsvAccess(gpu_manager).host_to_device_cb_idx_
    = CmiRegisterHandler((CmiHandler)hostToDeviceCallback);
  CsvAccess(gpu_manager).kernel_cb_idx_
    = CmiRegisterHandler((CmiHandler)kernelCallback);
  CsvAccess(gpu_manager).device_to_host_cb_idx_
    = CmiRegisterHandler((CmiHandler)deviceToHostCallback);
  CsvAccess(gpu_manager).light_cb_idx_
    = CmiRegisterHandler((CmiHandler)lightCallback);
#endif
}

#ifdef HAPI_CUDA_CALLBACK
// Callback function invoked by the CUDA runtime certain parts of GPU work are
// complete. It sends a converse message to the original PE to free the relevant
// device memory and invoke the user's callback. The reason for this method is
// that a thread created by the CUDA runtime does not have access to any of the
// CpvDeclare'd variables as it is not one of the threads created by the Charm++
// runtime.
static void CUDART_CB CUDACallback(cudaStream_t stream, cudaError_t status,
                                   void *data) {
#ifdef HAPI_NVTX_PROFILE
  NVTXTracer nvtx_range("CUDACallback", NVTXColor::Silver);
#endif

  if (status == cudaSuccess) {
    // send message to the original PE
    char *conv_msg = (char*)data;
    int dstRank = *((int *)(conv_msg + CmiMsgHeaderSizeBytes));
    CmiPushPE(dstRank, conv_msg);
  }
  else {
    CmiAbort("[HAPI] error before CUDACallback");
  }
}

enum CallbackStage {
  AfterHostToDevice,
  AfterKernel,
  AfterDeviceToHost
};

static void addCallback(hapiWorkRequest *wr, CallbackStage stage) {
  // create converse message to be delivered to this PE after CUDA callback
  char *conv_msg = (char *)CmiAlloc(CmiMsgHeaderSizeBytes + sizeof(int) +
                                  sizeof(hapiWorkRequest *)); // FIXME memory leak?
  *((int *)(conv_msg + CmiMsgHeaderSizeBytes)) = CmiMyRank();
  *((hapiWorkRequest **)(conv_msg + CmiMsgHeaderSizeBytes + sizeof(int))) = wr;

  int handlerIdx;
  switch (stage) {
    case AfterHostToDevice:
      handlerIdx = CsvAccess(gpu_manager).host_to_device_cb_idx_;
      break;
    case AfterKernel:
      handlerIdx = CsvAccess(gpu_manager).kernel_cb_idx_;
      break;
    case AfterDeviceToHost:
      handlerIdx = CsvAccess(gpu_manager).device_to_host_cb_idx_;
      break;
    default: // wrong type
      CmiFree(conv_msg);
      return;
  }
  CmiSetHandler(conv_msg, handlerIdx);

  // add callback into CUDA stream
  hapiCheck(cudaStreamAddCallback(wr->stream, CUDACallback, (void*)conv_msg, 0));
}
#endif // HAPI_CUDA_CALLBACK

/******************** DEPRECATED ********************/
// User calls this function to offload work to the GPU.
void hapiEnqueue(hapiWorkRequest* wr) {
#ifdef HAPI_NVTX_PROFILE
  NVTXTracer nvtx_range("enqueue", NVTXColor::Pomegranate);
#endif

#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).progress_lock_);
#endif

  // allocate device memory
  CsvAccess(gpu_manager).allocateBuffers(wr);

  // transfer data to device
  CsvAccess(gpu_manager).hostToDeviceTransfer(wr);

  // add host-to-device transfer callback
  if (wr->host_to_device_cb) {
    // while there is an ongoing workrequest, quiescence should not be detected
    // even if all PEs seem idle
    CmiAssert(hapiQdCreate);
    hapiQdCreate(1);

#ifdef HAPI_CUDA_CALLBACK
    addCallback(wr, AfterHostToDevice);
#else
    recordEvent(wr->stream, wr->host_to_device_cb, NULL);
#endif
  }

  // run kernel
  CsvAccess(gpu_manager).runKernel(wr);

  // add kernel callback
  if (wr->kernel_cb) {
    CmiAssert(hapiQdCreate);
    hapiQdCreate(1);

#ifdef HAPI_CUDA_CALLBACK
    addCallback(wr, AfterKernel);
#else
    recordEvent(wr->stream, wr->kernel_cb, NULL);
#endif
  }

  // transfer data to host
  CsvAccess(gpu_manager).deviceToHostTransfer(wr);

  // add device-to-host transfer callback
  CmiAssert(hapiQdCreate);
  hapiQdCreate(1);
#ifdef HAPI_CUDA_CALLBACK
  // always invoked to free memory
  addCallback(wr, AfterDeviceToHost);
#else
  if (wr->device_to_host_cb) {
    recordEvent(wr->stream, wr->device_to_host_cb, NULL, wr);
  }
  else {
    recordEvent(wr->stream, NULL, NULL, wr);
  }
#endif

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).progress_lock_);
#endif
}

/******************** DEPRECATED ********************/
// Creates a hapiWorkRequest object on the heap and returns it to the user.
hapiWorkRequest* hapiCreateWorkRequest() {
  return (new hapiWorkRequest);
}

hapiWorkRequest::hapiWorkRequest() :
    grid_dim(0), block_dim(0), shared_mem(0), host_to_device_cb(NULL),
    kernel_cb(NULL), device_to_host_cb(NULL), runKernel(NULL), state(0),
    user_data(NULL), free_user_data(false), free_host_to_device_cb(false),
    free_kernel_cb(false), free_device_to_host_cb(false)
  {
#ifdef HAPI_TRACE
    trace_name = "";
#endif
#ifdef HAPI_INSTRUMENT_WRS
    chare_index = -1;
#endif

#if CMK_SMP
    CmiLock(CsvAccess(gpu_manager).stream_lock_);
#endif

    // Create default per-PE streams if none exist
    if (CsvAccess(gpu_manager).getStream(0) == NULL) {
      CsvAccess(gpu_manager).createNStreams(CmiMyNodeSize());
    }

    stream = CsvAccess(gpu_manager).getStream(CmiMyRank() % CsvAccess(gpu_manager).getNStreams());

#if CMK_SMP
    CmiUnlock(CsvAccess(gpu_manager).stream_lock_);
#endif
  }


/******************** DEPRECATED ********************/
// Need to be updated with the Tracing API.
static inline void gpuEventStart(hapiWorkRequest* wr, int* index,
                                 WorkRequestStage event, ProfilingStage stage) {
#ifdef HAPI_TRACE
  gpuEventTimer* shared_gpu_events_ = CsvAccess(gpu_manager).gpu_events_;
  int shared_time_idx_ = CsvAccess(gpu_manager).time_idx_++;
  shared_gpu_events_[shared_time_idx_].cmi_start_time = CmiWallTimer();
  shared_gpu_events_[shared_time_idx_].event_type = event;
  shared_gpu_events_[shared_time_idx_].trace_name = wr->trace_name;
  *index = shared_time_idx_;
  shared_gpu_events_[shared_time_idx_].stage = stage;
#ifdef HAPI_DEBUG
  CmiPrintf("[HAPI] start event %d of WR %s, profiling stage %d\n",
         event, wr->trace_name, stage);
#endif
#endif // HAPI_TRACE
}

/******************** DEPRECATED ********************/
// Need to be updated with the Tracing API.
static inline void gpuEventEnd(int index) {
#ifdef HAPI_TRACE
  CsvAccess(gpu_manager).gpu_events_[index].cmi_end_time = CmiWallTimer();
  traceUserBracketEvent(CsvAccess(gpu_manager).gpu_events_[index].stage,
                        CsvAccess(gpu_manager).gpu_events_[index].cmi_start_time,
                        CsvAccess(gpu_manager).gpu_events_[index].cmi_end_time);
#ifdef HAPI_DEBUG
  Cmiprintf("[HAPI] end event %d of WR %s, profiling stage %d\n",
          CsvAccess(gpu_manager).gpu_events_[index].event_type,
          CsvAccess(gpu_manager).gpu_events_[index].trace_name,
          CsvAccess(gpu_manager).gpu_events_[index].stage);
#endif
#endif // HAPI_TRACE
}

static inline void hapiWorkRequestStartTime(hapiWorkRequest* wr) {
#ifdef HAPI_INSTRUMENT_WRS
  wr->phase_start_time = CmiWallTimer();
#endif
}

static inline void profileWorkRequestEvent(hapiWorkRequest* wr,
                                           WorkRequestStage event) {
#ifdef HAPI_INSTRUMENT_WRS
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).inst_lock_);
#endif

  if (CsvAccess(gpu_manager).init_instr_) {
    double tt = CmiWallTimer() - (wr->phase_start_time);
    int index = wr->chare_index;
    char type = wr->comp_type;
    char phase = wr->comp_phase;

    std::vector<hapiRequestTimeInfo> &vec = CsvAccess(gpu_manager).avg_times_[index][type];
    if (vec.size() <= phase) {
      vec.resize(phase+1);
    }
    switch (event) {
      case DataSetup:
        vec[phase].transfer_time += tt;
        break;
      case KernelExecution:
        vec[phase].kernel_time += tt;
        break;
      case DataCleanup:
        vec[phase].cleanup_time += tt;
        vec[phase].n++;
        break;
      default:
        CmiPrintf("[HAPI] invalid event during profileWorkRequestEvent\n");
    }
  }
  else {
    CmiPrintf("[HAPI] instrumentation not initialized!\n");
  }

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).inst_lock_);
#endif
#endif // HAPI_INSTRUMENT_WRS
}

// Create a pool with n_slots slots.
// There are n_buffers[i] buffers for each buffer size corresponding to entry i.
// TODO list the alignment/fragmentation issues with either of two allocation schemes:
// if single, large buffer is allocated for each subpool
// if multiple, smaller buffers are allocated for each subpool
static void createPool(int *n_buffers, int n_slots, std::vector<BufferPool> &pools){
  std::vector<size_t>& mempool_boundaries = CsvAccess(gpu_manager).mempool_boundaries_;

  // initialize pools
  pools.resize(n_slots);
  for (int i = 0; i < n_slots; i++) {
    pools[i].size = mempool_boundaries[i];
    pools[i].head = NULL;
  }

  int device;
  cudaDeviceProp device_prop;
  hapiCheck(cudaGetDevice(&device));
  hapiCheck(cudaGetDeviceProperties(&device_prop, device));

  // divide by # of PEs on physical node and multiply by # of PEs in logical node
  size_t available_memory = device_prop.totalGlobalMem /
                           CmiNumPesOnPhysicalNode(CmiPhysicalNodeID(CmiMyPe()))
                           * CmiMyNodeSize() * HAPI_MEMPOOL_SCALE;

  // pre-calculate memory per size
  int max_buffers = *std::max_element(n_buffers, n_buffers + n_slots);
  int n_buffers_to_allocate[n_slots];
  memset(n_buffers_to_allocate, 0, sizeof(n_buffers_to_allocate));
  size_t buf_size;
  while (available_memory >= mempool_boundaries[0] + sizeof(BufferPoolHeader)) {
    for (int i = 0; i < max_buffers; i++) {
      for (int j = n_slots - 1; j >= 0; j--) {
        buf_size = mempool_boundaries[j] + sizeof(BufferPoolHeader);
        if (i < n_buffers[j] && buf_size <= available_memory) {
          n_buffers_to_allocate[j]++;
          available_memory -= buf_size;
        }
      }
    }
  }

  // pin the host memory
  for (int i = 0; i < n_slots; i++) {
    buf_size = mempool_boundaries[i] + sizeof(BufferPoolHeader);
    int num_buffers = n_buffers_to_allocate[i];

    BufferPoolHeader* hd;
    BufferPoolHeader* previous = NULL;

    // pin host memory in a contiguous block for a slot
    void* pinned_chunk;
    hapiCheck(cudaMallocHost(&pinned_chunk, buf_size * num_buffers));

    // initialize header structs
    for (int j = num_buffers - 1; j >= 0; j--) {
      hd = reinterpret_cast<BufferPoolHeader*>(reinterpret_cast<unsigned char*>(pinned_chunk)
                                     + buf_size * j);
      hd->slot = i;
      hd->next = previous;
      previous = hd;
    }

    pools[i].head = previous;
#ifdef HAPI_MEMPOOL_DEBUG
    pools[i].num = num_buffers;
#endif
  }
}

static void releasePool(std::vector<BufferPool> &pools){
  for (int i = 0; i < pools.size(); i++) {
    BufferPoolHeader* hdr = pools[i].head;
    if (hdr != NULL) {
      hapiCheck(cudaFreeHost((void*)hdr));
    }
  }
  pools.clear();
}

static int findPool(size_t size){
  int boundary_array_len = CsvAccess(gpu_manager).mempool_boundaries_.size();
  if (size <= CsvAccess(gpu_manager).mempool_boundaries_[0]) {
    return 0;
  }
  else if (size > CsvAccess(gpu_manager).mempool_boundaries_[boundary_array_len-1]) {
    // create new slot
    CsvAccess(gpu_manager).mempool_boundaries_.push_back(size);

    BufferPool newpool;
    hapiCheck(cudaMallocHost((void**)&newpool.head, size + sizeof(BufferPoolHeader)));
    if (newpool.head == NULL) {
      CmiPrintf("[HAPI (%d)] findPool: failed to allocate newpool %d head, size %zu\n",
             CmiMyPe(), boundary_array_len, size);
      CmiAbort("[HAPI] failed newpool allocation");
    }
    newpool.size = size;
#ifdef HAPI_MEMPOOL_DEBUG
    newpool.num = 1;
#endif
    CsvAccess(gpu_manager).mempool_free_bufs_.push_back(newpool);

    BufferPoolHeader* hd = newpool.head;
    hd->next = NULL;
    hd->slot = boundary_array_len;

    return boundary_array_len;
  }
  for (int i = 0; i < CsvAccess(gpu_manager).mempool_boundaries_.size()-1; i++) {
    if (CsvAccess(gpu_manager).mempool_boundaries_[i] < size &&
        size <= CsvAccess(gpu_manager).mempool_boundaries_[i+1]) {
      return (i + 1);
    }
  }
  return -1;
}

static void* getBufferFromPool(int pool, size_t size){
  BufferPoolHeader* ret;

  if (pool < 0 || pool >= CsvAccess(gpu_manager).mempool_free_bufs_.size()) {
    CmiPrintf("[HAPI (%d)] getBufferFromPool, pool: %d, size: %zu invalid pool\n",
           CmiMyPe(), pool, size);
#ifdef HAPI_MEMPOOL_DEBUG
    CmiPrintf("[HAPI (%d)] num: %d\n", CmiMyPe(),
           CsvAccess(gpu_manager).mempool_free_bufs_[pool].num);
#endif
    CmiAbort("[HAPI] exiting after invalid pool");
  }
  else if (CsvAccess(gpu_manager).mempool_free_bufs_[pool].head == NULL) {
    BufferPoolHeader* hd;
    hapiCheck(cudaMallocHost((void**)&hd, sizeof(BufferPoolHeader) +
                             CsvAccess(gpu_manager).mempool_free_bufs_[pool].size));
#ifdef HAPI_MEMPOOL_DEBUG
    CmiPrintf("[HAPI (%d)] getBufferFromPool, pool: %d, size: %zu expand by 1\n",
           CmiMyPe(), pool, size);
#endif
    if (hd == NULL) {
      CmiAbort("[HAPI] exiting after NULL hd from pool");
    }
    hd->slot = pool;
    return (void*)(hd + 1);
  }
  else {
    ret = CsvAccess(gpu_manager).mempool_free_bufs_[pool].head;
    CsvAccess(gpu_manager).mempool_free_bufs_[pool].head = ret->next;
#ifdef HAPI_MEMPOOL_DEBUG
    ret->size = size;
    CsvAccess(gpu_manager).mempool_free_bufs_[pool].num--;
#endif
    return (void*)(ret + 1);
  }
  return NULL;
}

static void returnBufferToPool(int pool, BufferPoolHeader* hd) {
  hd->next = CsvAccess(gpu_manager).mempool_free_bufs_[pool].head;
  CsvAccess(gpu_manager).mempool_free_bufs_[pool].head = hd;
#ifdef HAPI_MEMPOOL_DEBUG
  CsvAccess(gpu_manager).mempool_free_bufs_[pool].num++;
#endif
}

void* hapiPoolMalloc(size_t size) {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).mempool_lock_);
#endif

  if (!CsvAccess(gpu_manager).mempool_initialized_) {
    // create pool of page-locked memory
    int sizes[HAPI_MEMPOOL_NUM_SLOTS];
          /*256*/ sizes[0]  =  4;
          /*512*/ sizes[1]  =  2;
         /*1024*/ sizes[2]  =  2;
         /*2048*/ sizes[3]  =  4;
         /*4096*/ sizes[4]  =  2;
         /*8192*/ sizes[5]  =  6;
        /*16384*/ sizes[6]  =  5;
        /*32768*/ sizes[7]  =  2;
        /*65536*/ sizes[8]  =  1;
       /*131072*/ sizes[9]  =  1;
       /*262144*/ sizes[10] =  1;
       /*524288*/ sizes[11] =  1;
      /*1048576*/ sizes[12] =  1;
      /*2097152*/ sizes[13] =  2;
      /*4194304*/ sizes[14] =  2;
      /*8388608*/ sizes[15] =  2;
     /*16777216*/ sizes[16] =  2;
     /*33554432*/ sizes[17] =  1;
     /*67108864*/ sizes[18] =  1;
    /*134217728*/ sizes[19] =  7;
    createPool(sizes, HAPI_MEMPOOL_NUM_SLOTS, CsvAccess(gpu_manager).mempool_free_bufs_);
    CsvAccess(gpu_manager).mempool_initialized_ = true;

#ifdef HAPI_MEMPOOL_DEBUG
    CmiPrintf("[HAPI (%d)] done creating buffer pool\n", CmiMyPe());
#endif
  }

  int pool = findPool(size);
  void* buf = getBufferFromPool(pool, size);

#ifdef HAPI_MEMPOOL_DEBUG
  CmiPrintf("[HAPI (%d)] hapiPoolMalloc size %zu pool %d left %d\n",
         CmiMyPe(), size, pool,
         CsvAccess(gpu_manager).mempool_free_bufs_[pool].num);
#endif

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).mempool_lock_);
#endif

  return buf;
}

void hapiPoolFree(void* ptr) {
  // check if mempool was initialized, just return if not
  if (!CsvAccess(gpu_manager).mempool_initialized_)
    return;

  BufferPoolHeader* hd = ((BufferPoolHeader*)ptr) - 1;
  int pool = hd->slot;

#ifdef HAPI_MEMPOOL_DEBUG
  size_t size = hd->size;
#endif

#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).mempool_lock_);
#endif

  returnBufferToPool(pool, hd);

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).mempool_lock_);
#endif

#ifdef HAPI_MEMPOOL_DEBUG
  CmiPrintf("[HAPI (%d)] hapiPoolFree size %zu pool %d left %d\n",
         CmiMyPe(), size, pool,
         CsvAccess(gpu_manager).mempool_free_bufs_[pool].num);
#endif
}

#ifdef HAPI_INSTRUMENT_WRS
void hapiInitInstrument(int n_chares, int n_types) {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).inst_lock_);
#endif

  if (!CsvAccess(gpu_manager).init_instr_) {
    CsvAccess(gpu_manager).avg_times_.resize(n_chares);
    for (int i = 0; i < n_chares; i++) {
      CsvAccess(gpu_manager).avg_times_[i].resize(n_types);
    }
    CsvAccess(gpu_manager).init_instr_ = true;
  }

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).inst_lock_);
#endif
}

hapiRequestTimeInfo* hapiQueryInstrument(int chare, char type, char phase) {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).inst_lock_);
#endif

  if (phase < CsvAccess(gpu_manager).avg_times_[chare][type].size()) {
    return &CsvAccess(gpu_manager).avg_times_[chare][type][phase];
  }
  else {
    return NULL;
  }

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).inst_lock_);
#endif
}

void hapiClearInstrument() {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).inst_lock_);
#endif

  for (int chare = 0; chare < CsvAccess(gpu_manager).avg_times_.size(); chare++) {
    for (char type = 0; type < CsvAccess(gpu_manager).avg_times_[chare].size(); type++) {
      CsvAccess(gpu_manager).avg_times_[chare][type].clear();
    }
    CsvAccess(gpu_manager).avg_times_[chare].clear();
  }
  CsvAccess(gpu_manager).avg_times_.clear();
  CsvAccess(gpu_manager).init_instr_ = false;

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).inst_lock_);
#endif
}
#endif // HAPI_INSTRUMENT_WRS

// Poll HAPI events stored in the PE's queue. Current strategy is to process
// all successive completed events in the queue starting from the front.
// TODO Maybe we should make one pass of all events in the queue instead,
// since there might be completed events later in the queue.
void hapiPollEvents() {
#ifndef HAPI_CUDA_CALLBACK
  std::queue<hapiEvent>& queue = CpvAccess(hapi_event_queue);
  while (!queue.empty()) {
    hapiEvent hev = queue.front();
    if (cudaEventQuery(hev.event) == cudaSuccess) {
      // invoke Charm++ callback if one was given
      if (hev.cb) {
        CmiAssert(hapiInvokeCallback);
        hapiInvokeCallback(hev.cb, hev.cb_msg);
      }

      // clean up hapiWorkRequest
      if (hev.wr) {
        hapiWorkRequestCleanup(hev.wr);
      }
      cudaEventDestroy(hev.event);
      queue.pop();
      CpvAccess(n_hapi_events)--;

      // inform QD that an event was processed
      CmiAssert(hapiQdProcess);
      hapiQdProcess(1);
    }
    else {
      // stop going through the queue once we encounter a non-successful event
      break;
    }
  }
#endif
}

int hapiCreateStreams() {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).stream_lock_);
#endif

  int ret = CsvAccess(gpu_manager).createStreams();

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).stream_lock_);
#endif

  return ret;
}

cudaStream_t hapiGetStream() {
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).stream_lock_);
#endif

  cudaStream_t ret = CsvAccess(gpu_manager).getNextStream();

#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).stream_lock_);
#endif

  return ret;
}

// Lightweight HAPI, to be invoked after data transfer or kernel execution.
void hapiAddCallback(cudaStream_t stream, void* cb, void* cb_msg) {
#ifndef HAPI_CUDA_CALLBACK
  // record CUDA event
  recordEvent(stream, cb, cb_msg);
#else
  /* FIXME works for now (faster too), but CmiAlloc might not be thread-safe
#if CMK_SMP
  CmiLock(CsvAccess(gpu_manager).queue_lock_);
#endif
*/

  // create converse message to be delivered to this PE after CUDA callback
  char* conv_msg = (char*)CmiAlloc(CmiMsgHeaderSizeBytes + sizeof(int) +
                                 sizeof(void*)); // FIXME memory leak?
  char* conv_msg_tmp = conv_msg + CmiMsgHeaderSizeBytes;
  *((int*)conv_msg_tmp) = CmiMyRank();
  conv_msg_tmp += sizeof(int);
  *((void**)conv_msg_tmp) = cb;
  CmiSetHandler(conv_msg, CsvAccess(gpu_manager).light_cb_idx_);

  // push into CUDA stream
  hapiCheck(cudaStreamAddCallback(stream, CUDACallback, (void*)conv_msg, 0));

  /*
#if CMK_SMP
  CmiUnlock(CsvAccess(gpu_manager).queue_lock_);
#endif
*/
#endif

  // while there is an ongoing workrequest, quiescence should not be detected
  // even if all PEs seem idle
  CmiAssert(hapiQdCreate);
  hapiQdCreate(1);
}

cudaError_t hapiMalloc(void** devPtr, size_t size) {
  return cudaMalloc(devPtr, size);
}

cudaError_t hapiFree(void* devPtr) {
  return cudaFree(devPtr);
}

cudaError_t hapiMallocHost(void** ptr, size_t size) {
  return cudaMallocHost(ptr, size);
}

cudaError_t hapiMallocHostPool(void** ptr, size_t size) {
  void* tmp_ptr = hapiPoolMalloc(size);
  if (tmp_ptr) {
    *ptr = tmp_ptr;
    return cudaSuccess;
  }
  else return cudaErrorMemoryAllocation;
}

cudaError_t hapiFreeHost(void* ptr) {
  return cudaFreeHost(ptr);
}

cudaError_t hapiFreeHostPool(void *ptr) {
  hapiPoolFree(ptr);
  return cudaSuccess;
}

cudaError_t hapiMemcpyAsync(void* dst, const void* src, size_t count, cudaMemcpyKind kind, cudaStream_t stream = 0) {
  return cudaMemcpyAsync(dst, src, count, kind, stream);
}

void hapiErrorDie(cudaError_t retCode, const char* code, const char* file, int line) {
  if (retCode != cudaSuccess) {
    fprintf(stderr, "Fatal CUDA Error [%d] %s at %s:%d\n", retCode, cudaGetErrorString(retCode), file, line);
    CmiAbort("Exit due to CUDA error");
  }
}
