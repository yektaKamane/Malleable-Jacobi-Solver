/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/**
 * \addtogroup CkPerf
*/
/*@{*/

#ifndef __TRACE_COMMON_H__
#define __TRACE_COMMON_H__

#define  CREATION           1
#define  BEGIN_PROCESSING   2
#define  END_PROCESSING     3
#define  ENQUEUE            4
#define  DEQUEUE            5
#define  BEGIN_COMPUTATION  6
#define  END_COMPUTATION    7
#define  BEGIN_INTERRUPT    8
#define  END_INTERRUPT      9
#define  MESSAGE_RECV       10
#define  BEGIN_TRACE        11
#define  END_TRACE          12
#define  USER_EVENT         13
#define  BEGIN_IDLE         14
#define  END_IDLE           15
#define  BEGIN_PACK         16
#define  END_PACK           17
#define  BEGIN_UNPACK       18
#define  END_UNPACK         19

#define  CREATION_MULTICAST 21

#define  USER_EVENT_PAIR    100

CkpvExtern(int, CtrLogBufSize);
CkpvExtern(char*, traceRoot);
CkpvExtern(double, traceInitTime);
CkpvExtern(double, traceInitCpuTime);


#if CMK_BLUEGENE_CHARM
#define  TRACE_TIMER   BgGetTime
#define  TRACE_CPUTIMER   BgGetTime
inline double TraceTimer() { return TRACE_TIMER(); }
inline double TraceTimer(double t) { return t; }
inline double TraceCpuTimer() { return TRACE_TIMER(); }
inline double TraceCpuTimer(double t) { return t; }
#else
#define  TRACE_TIMER   CmiWallTimer
#define  TRACE_CPUTIMER   CmiWallTimer
inline double TraceTimer() { return TRACE_TIMER() - CkpvAccess(traceInitTime); }
inline double TraceTimer(double t) { return t - CkpvAccess(traceInitTime); }
inline double TraceCpuTimer() { return TRACE_CPUTIMER() - CkpvAccess(traceInitCpuTime); }
inline double TraceCpuTimer(double t) { return t - CkpvAccess(traceInitCpuTime); }
#endif

/** Tracing-specific registered Charm entities: */
extern int _threadMsg, _threadChare, _threadEP;
extern int _packMsg, _packChare, _packEP;
extern int _unpackMsg, _unpackChare, _unpackEP;
extern int _dummyMsg, _dummyChare, _dummyEP;

/** Write out the common parts of the .sts file. */
extern void traceWriteSTS(FILE *stsfp,int nUserEvents);

#endif

/*@}*/
