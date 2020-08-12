#ifndef _MACHINE_RDMA_H_
#define _MACHINE_RDMA_H_

/*Function Pointer to Acknowledgement Handler*/
typedef void (*RdmaAckHandlerFn)(void *token);

/* Support for Nocopy Direct API */
typedef struct _cmi_common_rdma_info {
#if CMK_USE_CMA
  pid_t pid;
#elif defined _MSC_VER
  char empty;
#endif
} CmiCommonRdmaInfo_t;

/* Set the generic converse/LRTS information */
void CmiSetRdmaCommonInfo(void *info, const void *ptr, int size) {
#if CMK_USE_CMA
  CmiCommonRdmaInfo_t *cmmInfo = (CmiCommonRdmaInfo_t *)info;
  cmmInfo->pid = getpid();
#endif
}

int CmiGetRdmaCommonInfoSize() {
#if CMK_USE_CMA
  return sizeof(CmiCommonRdmaInfo_t);
#else
  return 0; // If CMK_USE_CMA is false, sizeof(CmiCommonRdmaInfo_t) is 1 (size of an empty structure in C++)
            // However, 0 is returned since CMK_COMMON_NOCOPY_DIRECT_BYTES is set to 0 when CMK_USE_CMA is false
            // because the offset (returned by CmiGetRdmaCommonInfoSize) should equal CMK_COMMON_NOCOPY_DIRECT_BYTES
#endif
}

#if CMK_USE_CMA
#include <sys/uio.h> // for struct iovec
extern int cma_works;
int readShmCma(pid_t, char*, char*, size_t);
int writeShmCma(pid_t, char *, char *, size_t);

// These methods are also used by the generic layer implementation of the Direct API
void CmiIssueRgetUsingCMA(
  const void* srcAddr,
  void *srcInfo,
  int srcPe,
  const void* destAddr,
  void *destInfo,
  int destPe,
  size_t size) {

  // get remote process id
  CmiCommonRdmaInfo_t *remoteCommInfo = (CmiCommonRdmaInfo_t *)srcInfo;
  pid_t pid = remoteCommInfo->pid;
  readShmCma(pid, (char *)destAddr, (char *)srcAddr, size);
}

void CmiIssueRputUsingCMA(
  const void* destAddr,
  void *destInfo,
  int destPe,
  const void* srcAddr,
  void *srcInfo,
  int srcPe,
  size_t size) {

  // get remote process id
  CmiCommonRdmaInfo_t *remoteCommInfo = (CmiCommonRdmaInfo_t *)destInfo;
  pid_t pid = remoteCommInfo->pid;
  writeShmCma(pid, (char *)srcAddr, (char *)destAddr, size);
}
#endif

#if CMK_ONESIDED_IMPL

// Function Pointer to the acknowledement handler function for the Direct API
RdmaAckHandlerFn ncpyDirectAckHandlerFn;

typedef struct _cmi_rdma_direct_ack {
  const void *srcAddr;
  int srcPe;
  const void *destAddr;
  int destPe;
  int ackSize;
} CmiRdmaDirectAck;

/* Support for Nocopy Direct API */
void LrtsSetRdmaBufferInfo(void *info, const void *ptr, int size, unsigned short int mode);
void LrtsSetRdmaNcpyAck(RdmaAckHandlerFn fn);
void LrtsIssueRget(NcpyOperationInfo *ncpyOpInfo);

void LrtsIssueRput(NcpyOperationInfo *ncpyOpInfo);

void LrtsDeregisterMem(const void *ptr, void *info, int pe, unsigned short int mode);

#if CMK_REG_REQUIRED
void LrtsInvokeRemoteDeregAckHandler(int pe, NcpyOperationInfo *ncpyOpInfo);
#endif

/* Set the machine specific information for a nocopy pointer */
void CmiSetRdmaBufferInfo(void *info, const void *ptr, int size, unsigned short int mode){
  LrtsSetRdmaBufferInfo(info, ptr, size, mode);
}

void CmiInvokeNcpyAck(void *ack) {
  ncpyDirectAckHandlerFn(ack);
}

/* Set the ack handler function used in the Direct API */
void CmiSetDirectNcpyAckHandler(RdmaAckHandlerFn fn){
  ncpyDirectAckHandlerFn = fn;
}

/* Perform an RDMA Get operation into the local destination address from the remote source address*/
void CmiIssueRget(NcpyOperationInfo *ncpyOpInfo) {
  // Use network RDMA for a PE on a remote host
  LrtsIssueRget(ncpyOpInfo);
}

/* Perform an RDMA Put operation into the remote destination address from the local source address */
void CmiIssueRput(NcpyOperationInfo *ncpyOpInfo) {
  // Use network RDMA for a PE on a remote host
  LrtsIssueRput(ncpyOpInfo);
}

/* De-register registered memory for pointer */
void CmiDeregisterMem(const void *ptr, void *info, int pe, unsigned short int mode){
  LrtsDeregisterMem(ptr, info, pe, mode);
}

#if CMK_REG_REQUIRED
void CmiInvokeRemoteDeregAckHandler(int pe, NcpyOperationInfo *ncpyOpInfo) {
  LrtsInvokeRemoteDeregAckHandler(pe, ncpyOpInfo);
}
#endif

#endif /*End of CMK_ONESIDED_IMPL */
#endif
