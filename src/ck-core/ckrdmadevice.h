#ifndef _CKRDMADEVICE_H_
#define _CKRDMADEVICE_H_

#include "ckcallback.h"
#include "conv-rdmadevice.h"

#if CMK_CUDA
#include <cuda_runtime.h>

#define CkNcpyModeDevice CmiNcpyModeDevice

struct CkDeviceBufferPost {
  // CUDA stream for device transfers
  cudaStream_t cuda_stream;

  // Use per-thread stream by default
  CkDeviceBufferPost() : cuda_stream(cudaStreamPerThread) {}
};

class CkDeviceBuffer : public CmiDeviceBuffer {
public:
  // Callback to be invoked on the sender/receiver
  CkCallback cb;

  CkDeviceBuffer() : CmiDeviceBuffer() {}

  explicit CkDeviceBuffer(const void* ptr_) : CmiDeviceBuffer(ptr_, 0) {
    cb = CkCallback(CkCallback::ignore);
  }

  explicit CkDeviceBuffer(const void* ptr_, const CkCallback& cb_) : CmiDeviceBuffer(ptr_, 0) {
    cb = cb_;
  }

  explicit CkDeviceBuffer(const void* ptr_, cudaStream_t cuda_stream_) : CmiDeviceBuffer(ptr_, 0) {
    cb = CkCallback(CkCallback::ignore);
    cuda_stream = cuda_stream_;
  }

  explicit CkDeviceBuffer(const void* ptr_, const CkCallback& cb_, cudaStream_t cuda_stream_) : CmiDeviceBuffer(ptr_, 0) {
    cb = cb_;
    cuda_stream = cuda_stream_;
  }


  explicit CkDeviceBuffer(const void* ptr_, size_t cnt_) : CmiDeviceBuffer(ptr_, cnt_) {
    cb = CkCallback(CkCallback::ignore);
  }

  explicit CkDeviceBuffer(const void* ptr_, size_t cnt_, const CkCallback& cb_) : CmiDeviceBuffer(ptr_, cnt_) {
    cb = cb_;
  }

  explicit CkDeviceBuffer(const void* ptr_, size_t cnt_, cudaStream_t cuda_stream_) : CmiDeviceBuffer(ptr_, cnt_) {
    cb = CkCallback(CkCallback::ignore);
    cuda_stream = cuda_stream_;
  }

  explicit CkDeviceBuffer(const void* ptr_, size_t cnt_, const CkCallback& cb_, cudaStream_t cuda_stream_) : CmiDeviceBuffer(ptr_, cnt_) {
    cb = cb_;
    cuda_stream = cuda_stream_;
  }

  void pup(PUP::er &p) {
    CmiDeviceBuffer::pup(p);
    p|cb;
  }

  friend bool CkRdmaDeviceIssueRgets(envelope *env, int numops, void **arrPtrs, int *arrSizes, CkDeviceBufferPost *postStructs);
};

bool CkRdmaDeviceIssueRgets(envelope *env, int numops, void **arrPtrs, int *arrSizes, CkDeviceBufferPost *postStructs);
void CkRdmaDeviceOnSender(int dest_pe, int numops, CkDeviceBuffer** buffers);

#endif // CMK_CUDA

#endif // _CKRDMADEVICE_H_
