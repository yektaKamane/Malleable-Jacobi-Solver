/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/
/**
\file
\brief Charm Kernel--the groups and chares level of Charm++.
*/
#ifndef CHARM_H
#define CHARM_H

#include "converse.h"
#include <sys/types.h> /* for size_t */

#ifdef __cplusplus
#include "pup.h"
extern "C" {
#endif

/******************************************************************************
 *
 * Converse Concepts, renamed to CK
 *
 *****************************************************************************/

/** Queueing types, for use with CkSetQueueing: */
#define CK_QUEUEING_FIFO   CQS_QUEUEING_FIFO
#define CK_QUEUEING_LIFO   CQS_QUEUEING_LIFO
#define CK_QUEUEING_IFIFO  CQS_QUEUEING_IFIFO
#define CK_QUEUEING_ILIFO  CQS_QUEUEING_ILIFO
#define CK_QUEUEING_BFIFO  CQS_QUEUEING_BFIFO
#define CK_QUEUEING_BLIFO  CQS_QUEUEING_BLIFO

#define CkTimer  	CmiTimer
#define CkWallTimer  	CmiWallTimer
#define CkCpuTimer  	CmiCpuTimer

#define CkMyPe		CmiMyPe
#define CkMyRank	CmiMyRank
#define CkMyNode	CmiMyNode
#define CkNumPes	CmiNumPes
#define CkNumNodes	CmiNumNodes
#define CkNodeFirst	CmiNodeFirst
#define CkNodeSize	CmiNodeSize
#define CkMyNodeSize	CmiMyNodeSize
#define CkNodeOf	CmiNodeOf
#define CkRankOf	CmiRankOf

#define CkPrintf                CmiPrintf
#define CkScanf                 CmiScanf
#define CkError                 CmiError
#define CkAbort                 CmiAbort
#define CkAssert                CmiAssert
extern void  CkExit(void);
extern char **CkGetArgv(void);
extern int  CkGetArgc(void);

/******************************************************************************
 *
 * Miscellaneous Constants
 *
 *****************************************************************************/

#define CK_PE_ALL        CLD_BROADCAST_ALL
#define CK_PE_ALL_BUT_ME CLD_BROADCAST
#define CK_PE_ANY        CLD_ANYWHERE

/******************************************************************************
 *
 * Message Allocation Calls
 *
 *****************************************************************************/

extern void* CkAllocSysMsg(void);
extern void  CkFreeSysMsg(void *msg);
extern void* CkAllocMsg(int msgIdx, int msgBytes, int prioBits);
extern void* CkAllocBuffer(void *msg, int bufsize);
extern void  CkFreeMsg(void *msg);
extern void* CkCopyMsg(void **pMsg);
extern void  CkSetQueueing(void *msg, int strategy);
extern void* CkPriorityPtr(void *msg);

/*********************************************************/
/**
\addtogroup CkRegister
\brief Charm Registration--keeps track of the possible chare and method types.

These are implemented in register.C.
*/
/*@{*/
/** Message pack function: convert a message into a buffer. */
typedef void* (*CkPackFnPtr)(void *msg);
/** Message unpack function: convert a buffer into a message. */
typedef void* (*CkUnpackFnPtr)(void *buf);

/** Register this message name, with this basic size and pack and unpack functions. */
extern int CkRegisterMsg(const char *name, CkPackFnPtr pack, 
                       CkUnpackFnPtr unpack, size_t size);

/** This entry point flag indicates the method does not keep the passed-in message. */
#define CK_EP_NOKEEP       (1<<2) 
#define CK_EP_INTRINSIC    (1<<3) 
#define CK_EP_TRACEDIABLE  (1<<4) 

/** A "call function" to invoke a method on an object. See EntryInfo */
typedef void  (*CkCallFnPtr) (void *msg, void *obj);
/** Register this entry point, with this call function and flags.
    Returns the entry point's index in the _entryTable. */
extern int CkRegisterEp(const char *name, CkCallFnPtr call, int msgIdx, 
                        int chareIdx,int ck_ep_flags);

/** Register this type of chare (group, or array), with this size.
    Returns the Chare's index in the _chareTable. */
extern int CkRegisterChare(const char *name, int dataSz);
/** Register this chare as a mainchare, with this entry point as its constructor.*/
extern int CkRegisterMainChare(int chareIndex, int epIndex);
/** Register a default constructor for this chare.*/
extern void CkRegisterDefaultCtor(int chareIndex, int ctorEpIndex);
/** Register a migration constructor for this chare.*/
extern void CkRegisterMigCtor(int chareIndex, int ctorEpIndex);
/** Indicate whether this group is an IrrGroup. */
extern void CkRegisterGroupIrr(int chareIndex,int isIrr);
/** Register the chare baseIdx as a base class of the chare derivedIdx. */
extern void CkRegisterBase(int derivedIdx, int baseIdx);

/** This function pup's a global variable.*/
typedef void (*CkPupReadonlyFnPtr)(void *pup_er);
/** Register this readonly global variable.*/
extern void CkRegisterReadonly(const char *name,const char *type,
	int size, void *ptr,CkPupReadonlyFnPtr pup_fn);
/** Register this readonly message.*/
extern void CkRegisterReadonlyMsg(const char *name,const char *type,
	void** pMsg);

/** A "marshall unpack" function: pups out parameters and calls a method. */
typedef int (*CkMarshallUnpackFn)(char *marshall_buf,void *object);
/** Register this marshall unpack function with this entry point.*/
extern void CkRegisterMarshallUnpackFn(int epIndex,CkMarshallUnpackFn m);
/** Lookup the marshall unpack function, if any, for this entry point.*/
extern CkMarshallUnpackFn CkLookupMarshallUnpackFn(int epIndex);

#ifdef __cplusplus
/** A "message pup" function: pups message data for debugger display. */
typedef void (*CkMessagePupFn)(PUP::er &p,void *userMessage);
/** Register this message pup function with this entry point.*/
extern void CkRegisterMessagePupFn(int epIndex,CkMessagePupFn m);
#endif
/*@}*/

/*********************************************************/
/**
\addtogroup Ck
\brief Charm Kernel--the groups and chares level of Charm++.

These routines are implemented in ck.C.
*/
/*@{*/

typedef struct {
  int   onPE;
  void* objPtr;
} CkChareID;

typedef struct _ckGroupID{
  int idx;		/* pe(processor number) is removed from the structure */
#ifdef __cplusplus
  void pup(PUP::er &p) {  p|idx; }
  int isZero(void) const { return (idx==0); }
  void setZero(void) { idx=0; }
  int operator==(const struct _ckGroupID& gid) const {
    return (gid.idx==idx);
  }
#endif
} CkGroupID;

typedef CkGroupID CkNodeGroupID;

/******************************************************************************
 *
 * Object Creation Calls
 *
 *****************************************************************************/
#ifdef __cplusplus
class envelope;
#else
typedef struct envelope envelope;
#endif
extern void CkCreateChare(int chareIdx, int constructorIdx, void *msg,
                          CkChareID *vid, int destPE);
extern CkGroupID CkCreateGroup(int chareIdx, int constructorIdx, void *msg);
extern CkGroupID CkCreateNodeGroup(int chareIdx, int constructorIdx, void *msg);
extern void CkCreateLocalGroup(CkGroupID groupID, int constructorIdx, envelope *env);
extern void CkCreateLocalNodeGroup(CkGroupID groupID, int constructorIdx, envelope *env);

/******************************************************************************
 *
 * Asynchronous Remote Method Invocation Calls
 *
 *****************************************************************************/

extern void CkSendMsg(int entryIndex, void *msg, const CkChareID *chare);
extern void CkSendMsgBranch(int eIdx, void *msg, int destPE, CkGroupID gID);
extern void CkSendMsgInline(int entryIndex, void *msg, const CkChareID *chare);
extern void CkSendMsgBranchInline(int eIdx, void *msg, int destPE, CkGroupID gID);
extern void CkSendMsgBranchMulti(int eIdx, void *msg, int npes, int *pes, 
                                 CkGroupID gID);
extern void CkSendMsgNodeBranch(int eIdx, void *msg, int destNode, 
                                CkGroupID gID);
extern void CkSendMsgNodeBranchInline(int eIdx, void *msg, int destNode, 
                                CkGroupID gID);
extern void CkBroadcastMsgBranch(int eIdx, void *msg, CkGroupID gID);
extern void CkBroadcastMsgNodeBranch(int eIdx, void *msg, CkGroupID gID);

extern int  CkChareMsgPrep(int eIdx, void *msg,const CkChareID *pCid);
extern void CkGroupMsgPrep(int eIdx, void *msg, CkGroupID gID);
extern void CkNodeGroupMsgPrep(int eIdx, void *msg, CkGroupID gID);

extern void CkSetRefNum(void *msg, int ref);
extern int  CkGetRefNum(void *msg);
extern int  CkGetSrcPe(void *msg);
extern int  CkGetSrcNode(void *msg);

extern void CkDeliverMessageFree(int epIdx,void *msg,void *object);
extern void CkDeliverMessageReadonly(int epIdx,const void *msg,void *object);

extern void *CkLocalBranch(CkGroupID gID);
extern void *CkLocalNodeBranch(CkGroupID gID);
extern void *CkLocalChare(const CkChareID *chare);

extern void CkArrayManagerInsert(int onPe,void *msg,CkGroupID aID);
extern void CkArrayManagerDeliver(int onPe,void *msg);

/*@}*/



/*********************************************************/
/**
\addtogroup CkFutures
\brief Futures--ways to block Converse threads on remote events.

These routines are implemented in ckfutures.C.
*/
/*@{*/
typedef int CkFutureID;
extern void* CkRemoteCall(int eIdx, void *msg,const CkChareID *chare);
extern void* CkRemoteBranchCall(int eIdx, void *msg, CkGroupID gID, int pe);
extern void* CkRemoteNodeBranchCall(int eIdx, void *msg, CkGroupID gID, int node);
extern CkFutureID CkRemoteCallAsync(int eIdx, void *msg, const CkChareID *chare);
extern CkFutureID CkRemoteBranchCallAsync(int eIdx, void *msg, CkGroupID gID, 
                                          int pe);
extern CkFutureID CkRemoteNodeBranchCallAsync(int eIdx, void *msg, 
                                              CkGroupID gID, int node);

extern void* CkWaitFuture(CkFutureID futNum);
extern void CkWaitVoidFuture(CkFutureID futNum);
extern void CkReleaseFuture(CkFutureID futNum);
extern int CkProbeFuture(CkFutureID futNum);
extern void  CkSendToFuture(CkFutureID futNum, void *msg, int pe);
extern CkFutureID CkCreateAttachedFuture(void *msg);
extern void *CkWaitReleaseFuture(CkFutureID futNum);

/******************************************************************************
 *
 * Semaphore calls
 *
 *****************************************************************************/

typedef struct _ckSemaID {
  int pe;
  int idx;
#ifdef __cplusplus
  public:
    void pup(PUP::er &p) { p(pe); p(idx); }
#endif
} CkSemaID;

extern CkSemaID CkSemaCreate(void);
extern void *CkSemaWait(CkSemaID id);
extern void CkSemaWaitN(CkSemaID id, int n, void *marray[]);
extern void CkSemaSignal(CkSemaID id, void *m);
extern void CkSemaDestroy(CkSemaID id);
/*@}*/


/******************************************************************************
 *
 * Quiescence Calls
 *
 *****************************************************************************/
/**
\addtogroup CkQD
\brief Quiescence Detection--a way to tell when nothing is happening.

These routines are implemented in qd.C and waitqd.C.
*/
/*@{*/

/** When quiescence occurs, send a message to this entry point of this Chare. */
extern void CkStartQD(int eIdx,const CkChareID *chare);
/** Block until quiescence occurs. */
extern void CkWaitQD(void);
/*@}*/

/******************************************************************************
 *
 * Miscellaneous Calls
 *
 *****************************************************************************/

extern int CkMessageToEpIdx(void *msg);
extern void CkPrintEntryMethod(int epIdx);
extern void CkPrintChareName(int chareIdx);
extern void CkSummary_MarkEvent(int);
extern void CkSummary_StartPhase(int);
extern int CkDisableTracing(int epIdx);
extern void CkEnableTracing(int epIdx);

#ifdef __cplusplus
}
#endif
#endif
