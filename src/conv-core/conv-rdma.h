#ifndef _CONV_RDMA_H
#define _CONV_RDMA_H

typedef void (*RdmaSingleAckCallerFn)(void *cbPtr, int pe, const void *ptr);
typedef void (*RdmaAckCallerFn)(void *token);

void *CmiSetRdmaAck(RdmaAckCallerFn fn, void *token);
void CmiSetRdmaInfo(void *dest, int destPE, int numOps);
void CmiSetRdmaOpInfo(void *dest, const void *ptr, int size, void *ack, int destPE);
int CmiGetRdmaOpInfoSize(void);
int CmiGetRdmaGenInfoSize(void);

int CmiGetRdmaInfoSize(int numOps);
void CmiSetRdmaRecvInfo(void *dest, int numOps, void *msg, void *rdmaInfo, int msgSize);
void CmiSetRdmaRecvOpInfo(void *dest, void *buffer, void *src_ref, int size, int opIndex, void *rdmaInfo);
int CmiGetRdmaOpRecvInfoSize(void);
int CmiGetRdmaGenRecvInfoSize(void);
int CmiGetRdmaRecvInfoSize(int numOps);

void CmiIssueRgets(void *recv, int pe);

/* Support for Direct API */
void CmiSetRdmaCommonInfo(void *info, const void *ptr, int size);
int CmiGetRdmaCommonInfoSize(void);

void CmiSetRdmaSrcInfo(void *info, const void *ptr, int size, unsigned short int mode);
void CmiSetRdmaDestInfo(void *info, const void *ptr, int size, unsigned short int mode);
void CmiSetRdmaNcpyAck(RdmaSingleAckCallerFn fn);

/* CmiIssueRget initiates an RDMA read operation, transferring 'size' bytes of data from the address space of 'srcPe' to local address, 'destAddr'.
 * When the runtime invokes srcAck on the source (target), it indicates safety to overwrite or free the srcAddr buffer.
 * When the runtime invokes destAck on the destination (initiator), it indicates that the data has been successfully received in the
 * destAddr buffer.
 */
void CmiIssueRget(
  const void* srcAddr,
  void *srcInfo,
  void *srcAck,
  int srcAckSize,
  int srcPe,
  unsigned short int *srcMode,
  const void* destAddr,
  void *destInfo,
  void *destAck,
  int destAckSize,
  int destPe,
  unsigned short int *destMode,
  int size);

/* CmiIssueRput initiates an RDMA write operation, transferring 'size' bytes of data from the local address, 'srcAddr' to the address space of 'destPe'.
 * When the runtime invokes srcAck on the source (initiator), it indicates safety to overwrite or free the srcAddr buffer.
 * When the runtime invokes destAck on the destination (target), it indicates that the data has been successfully received in the
 * destAddr buffer.
 */

void CmiIssueRput(
  const void* destAddr,
  void *destInfo,
  void *destAck,
  int destAckSize,
  int destPe,
  unsigned short int *destMode,
  const void* srcAddr,
  void *srcInfo,
  void *srcAck,
  int srcAckSize,
  int srcPe,
  unsigned short int *srcMode,
  int size);

void CmiReleaseSourceResource(const void *ptr, void *info, int pe, unsigned short int mode);
void CmiReleaseDestinationResource(const void *ptr, void *info, int pe, unsigned short int mode);

#if CMK_USE_CMA
void CmiIssueRgetUsingCMA(
  const void* srcAddr,
  void *srcInfo,
  int srcPe,
  const void* destAddr,
  void *destInfo,
  int destPe,
  int size);

void CmiIssueRputUsingCMA(
  const void* destAddr,
  void *destInfo,
  int destPe,
  const void* srcAddr,
  void *srcInfo,
  int srcPe,
  int size);
#endif

// Allocation from pool
void *CmiRdmaAlloc(int size);

#if !CMK_ONESIDED_DIRECT_IMPL
// Function declaration used for the generic implementation of the Nocopy Direct API
void CmiOnesidedDirectInit(void);
#endif

// Macros required to keep the Nocopy Direct API functional on non-LRTS layers
#if !CMK_USE_LRTS
#define CMK_BUFFER_REG                 0
#define CMK_BUFFER_UNREG               1
#define CMK_BUFFER_PREREG              2

#define CMK_COMMON_NOCOPY_DIRECT_BYTES 0
#endif

#endif
