#ifndef _MIDDLE_BLUE_H_
#define _MIDDLE_BLUE_H_

#include <memory.h>
#include "converse.h"
#include "blue.h"

#undef CkMyPe
#undef CkNumPes
#undef CkMyRank
#undef CkMyNode
#undef CkNumNodes
#undef CkMyNodeSize

#undef CmiSyncSend
#undef CmiSyncSendAndFree
#undef CmiSyncBroadcast
#undef CmiSyncBroadcastAndFree
#undef CmiSyncBroadcastAll
#undef CmiSyncBroadcastAllAndFree

#undef CmiSyncNodeSend
#undef CmiSyncNodeSendAndFree
#undef CmiSyncNodeBroadcast
#undef CmiSyncNodeBroadcastAndFree
#undef CmiSyncNodeBroadcastAll
#undef CmiSyncNodeBroadcastAllAndFree


#define CkVTimer   BgGetTime
#define CkElapse   BgElapse

#define CkRegisterHandler(x)        BgRegisterHandler((BgHandler)(x))
#define CkNumberHandler(n, x)       BgNumberHandler(n, (BgHandler)(x))
#define CkNumberHandlerEx(n, x, p)  BgNumberHandlerEx(n, (BgHandlerEx)(x), p)

#define ConverseExit             BgCharmExit

/**
  This version Blue Gene Charm++ use a whole Blue Gene node as 
  a Charm PE.
*/
#if CMK_BLUEGENE_NODE

#define CkpvDeclare 	BnvDeclare
#define CkpvExtern 	BnvExtern
#define CkpvStaticDeclare  BnvStaticDeclare
#define CkpvInitialize 	BnvInitialize
#define CkpvAccess	BnvAccess
#define CkpvAccessOther	BnvAccessOther

namespace BGConverse {

inline int CkMyPe() { return BgMyNode(); }
inline int CkNumPes() { int x,y,z; BgGetSize(&x, &y, &z); return (x*y*z); }
inline int CkMyRank() { return 0; }
inline int BgNodeRank() { return BgMyRank(); }
inline int CkMyNodeSize() { return 1; }

static inline void CmiSyncSend(int pe, int nb, char *m) 
{
  int x,y,z;
  char *dupm = (char *)CmiAlloc(nb);

//CmiPrintf("[%d] CmiSyncSend handle:%d\n", CkMyPe(), CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgGetXYZ(pe, &x, &y, &z);
  BgSendPacket(x,y,z, ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncSendAndFree(int pe, int nb, char *m)
{
  int x,y,z;
//CmiPrintf("[%d] CmiSyncSendAndFree handle:%d\n", CkMyPe(), CmiGetHandler(m));
  BgGetXYZ(pe, &x, &y, &z);
  BgSendPacket(x,y,z, ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncBroadcast(int nb, char *m)
{
  char *dupm = (char *)CmiAlloc(nb);
//CmiPrintf("[%d] CmiSyncBroadcast handle:%d\n", CkMyPe(), CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgBroadcastPacketExcept(CkMyPe(), ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncBroadcastAndFree(int nb, char *m)
{
//CmiPrintf("CmiSyncBroadcastAndFree handle:%d\n", CmiGetHandler(m));
  BgBroadcastPacketExcept(CkMyPe(), ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncBroadcastAll(int nb, char *m)
{
  char *dupm = (char *)CmiAlloc(nb);
//CmiPrintf("CmiSyncBroadcastAll: handle:%d\n", CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgBroadcastAllPacket(CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncBroadcastAllAndFree(int nb, char *m)
{
//CmiPrintf("CmiSyncBroadcastAllAndFree: handle:%d\n", CmiGetHandler(m));
  /* broadcast to all nodes */
  BgBroadcastAllPacket(CmiGetHandler(m), LARGE_WORK, nb, m);
}

}  /* end of namespace */



#else
/**
  This version of Blue Gene Charm++ use a Blue Gene thread as 
  a Charm PE.
*/

#define CkpvDeclare 	   BpvDeclare
#define CkpvExtern 	   BpvExtern
#define CkpvStaticDeclare  BpvStaticDeclare
#define CkpvInitialize 	   BpvInitialize
#define CkpvAccess	   BpvAccess
#define CkpvAccessOther	   BpvAccessOther

#define CksvDeclare 	   BnvDeclare
#define CksvExtern 	   BnvExtern
#define CksvStaticDeclare  BnvStaticDeclare
#define CksvInitialize 	   BnvInitialize
#define CksvAccess	   BnvAccess

namespace BGConverse {

static inline int CkMyPe() { return BgGetGlobalWorkerThreadID(); }
static inline int CkNumPes() { return BgNumNodes()*BgGetNumWorkThread(); }
static inline int CkMyRank() { return BgGetThreadID(); }
static inline int BgNodeRank() { return BgMyRank()*BgGetNumWorkThread()+BgGetThreadID(); }
static inline int CkMyNode() { return BgMyNode(); }
static inline int CkNumNodes() { return BgNumNodes(); }
static inline int CkMyNodeSize() { return BgGetNumWorkThread(); }

static inline void CmiSyncSend(int pe, int nb, char *m) 
{
  int x,y,z,t;
  char *dupm = (char *)CmiAlloc(nb);

//CmiPrintf("[%d] CmiSyncSend handle:%d\n", CkMyPe(), CmiGetHandler(m));
  memcpy(dupm, m, nb);
  t = pe%BgGetNumWorkThread();
  pe = pe/BgGetNumWorkThread();
  BgGetXYZ(pe, &x, &y, &z);
  BgSendPacket(x,y,z, t, CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncSendAndFree(int pe, int nb, char *m)
{
  int x,y,z,t;
//CmiPrintf("[%d] CmiSyncSendAndFree handle:%d\n", CkMyPe(), CmiGetHandler(m));
  t = pe%BgGetNumWorkThread();
  pe = pe/BgGetNumWorkThread();
  BgGetXYZ(pe, &x, &y, &z);
  BgSendPacket(x,y,z, t, CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncBroadcast(int nb, char *m)
{
  char *dupm = (char *)CmiAlloc(nb);
//CmiPrintf("[%d] CmiSyncBroadcast handle:%d\n", CkMyPe(), CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgThreadBroadcastPacketExcept(BgMyNode(), BgGetThreadID(), CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncBroadcastAndFree(int nb, char *m)
{
//CmiPrintf("CmiSyncBroadcastAndFree handle:%d node:%d tid:%d\n", CmiGetHandler(m), BgMyNode(), BgGetThreadID());
  BgThreadBroadcastPacketExcept(BgMyNode(), BgGetThreadID(), CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncBroadcastAll(int nb, char *m)
{
  char *dupm = (char *)CmiAlloc(nb);
//CmiPrintf("CmiSyncBroadcastAll: handle:%d\n", CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgThreadBroadcastAllPacket(CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncBroadcastAllAndFree(int nb, char *m)
{
//CmiPrintf("CmiSyncBroadcastAllAndFree: handle:%d\n", CmiGetHandler(m));
  /* broadcast to all nodes */
  BgThreadBroadcastAllPacket(CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncNodeSend(int node, int nb, char *m)
{
  int x,y,z,t;
  char *dupm = (char *)CmiAlloc(nb);

//CmiPrintf("[%d] CmiSyncNodeSend handle:%d\n", CkMyPe(), CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgGetXYZ(node, &x, &y, &z);
  BgSendPacket(x,y,z, ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncNodeSendAndFree(int node, int nb, char *m)
{
  int x,y,z,t;
//CmiPrintf("[%d] CmiSyncNodeSendAndFree handle:%d\n", CkMyPe(), CmiGetHandler(m));
  BgGetXYZ(node, &x, &y, &z);
  BgSendPacket(x,y,z, ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncNodeBroadcast(int nb, char *m)
{
  char *dupm = (char *)CmiAlloc(nb);
//CmiPrintf("[%d] CmiSyncBroadcast handle:%d\n", CkMyPe(), CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgBroadcastPacketExcept(CkMyNode(), ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncNodeBroadcastAndFree(int nb, char *m)
{
//CmiPrintf("CmiSyncBroadcastAndFree handle:%d node:%d tid:%d\n", CmiGetHandler(m), BgMyNode(), BgGetThreadID());
  BgBroadcastPacketExcept(CkMyNode(), ANYTHREAD, CmiGetHandler(m), LARGE_WORK, nb, m);
}

static inline void CmiSyncNodeBroadcastAll(int nb, char *m)
{
  char *dupm = (char *)CmiAlloc(nb);
//CmiPrintf("CmiSyncBroadcastAll: handle:%d\n", CmiGetHandler(m));
  memcpy(dupm, m, nb);
  BgBroadcastAllPacket(CmiGetHandler(m), LARGE_WORK, nb, dupm);
}

static inline void CmiSyncNodeBroadcastAllAndFree(int nb, char *m)
{
//CmiPrintf("CmiSyncBroadcastAllAndFree: handle:%d\n", CmiGetHandler(m));
  /* broadcast to all nodes */
  BgBroadcastAllPacket(CmiGetHandler(m), LARGE_WORK, nb, m);
}

}  /* end of namespace */

#endif


/** common functions for two versions */
namespace BGConverse {

static inline void BgCharmExit()
{
//  traceCharmClose();
  if (CkMyPe() == 0)  BgShutdown();
}

}


#endif
