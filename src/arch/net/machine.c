/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** @file
 * Basic NET implementation of Converse machine layer
 * @ingroup NET
 */

/** @defgroup NET
 * NET implementation of machine layer, ethernet in particular
 * @ingroup Machine
 *
 * THE DATAGRAM STREAM
 *
 * Messages are sent using UDP datagrams.  The sender allocates a
 * struct for each datagram to be sent.  These structs stick around
 * until slightly after the datagram is acknowledged.
 *
 * Datagrams are transmitted node-to-node (as opposed to pe-to-pe).
 * Each node has an OtherNode struct for every other node in the
 * system.  The OtherNode struct contains:
 *
 *   send_queue   (all datagram-structs not yet transmitted)
 *   send_window  (all datagram-structs transmitted but not ack'd)
 *
 * When an acknowledgement comes in, all packets in the send-window
 * are either marked as acknowledged or pushed back into the send
 * queue for retransmission.
 *
 * THE OUTGOING MESSAGE
 *
 * When you send or broadcast a message, the first thing the system
 * does is system creates an OutgoingMsg struct to represent the
 * operation.  The OutgoingMsg contains a very direct expression
 * of what you want to do:
 *
 * OutgoingMsg:
 *
 *   size      --- size of message in bytes
 *   data      --- pointer to the buffer containing the message
 *   src       --- processor which sent the message
 *   dst       --- destination processor (-1=broadcast, -2=broadcast all)
 *   freemode  --- see below.
 *   refcount  --- see below.
 *
 * The OutgoingMsg is kept around until the transmission is done, then
 * it is garbage collected --- the refcount and freemode fields are
 * to assist garbage collection.
 *
 * The freemode indicates which kind of buffer-management policy was
 * used (sync, async, or freeing).  The sync policy is handled
 * superficially by immediately converting sync sends into freeing
 * sends.  Thus, the freemode can either be 'A' (async) or 'F'
 * (freeing).  If the freemode is 'F', then garbage collection
 * involves freeing the data and the OutgoingMsg structure itself.  If
 * the freemode is 'A', then the only cleanup is to change the
 * freemode to 'X', a condition which is then detectable by
 * CmiAsyncMsgSent.  In this case, the actual freeing of the
 * OutgoingMsg is done by CmiReleaseCommHandle.
 *
 * When the transmission is initiated, the system computes how many
 * datagrams need to be sent, total.  This number is stored in the
 * refcount field.  Each time a datagram is delivered, the refcount
 * is decremented, when it reaches zero, cleanup is performed.  There
 * are two exceptions to this rule.  Exception 1: if the OutgoingMsg
 * is a send (not a broadcast) and can be performed with shared
 * memory, the entire datagram system is bypassed, the message is
 * simply delivered and freed, not using the refcount mechanism at
 * all.  Exception 2: If the message is a broadcast, then part of the
 * broadcast that can be done via shared memory is performed prior to
 * initiating the datagram/refcount system.
 *
 * DATAGRAM FORMATS AND MESSAGE FORMATS
 *
 * Datagrams have this format:
 *
 *   srcpe   (16 bits) --- source processor number.
 *   magic   ( 8 bits) --- magic number to make sure DG is good.
 *   dstrank ( 8 bits) --- destination processor rank.
 *   seqno   (32 bits) --- packet sequence number.
 *   data    (XX byte) --- user data.
 *
 * The only reason the srcpe is in there is because the receiver needs
 * to know which receive window to use.  The dstrank field is needed
 * because transmission is node-to-node.  Once the message is
 * assembled by the node, it must be delivered to the appropriate PE.
 * The dstrank field is used to encode certain special-case scenarios.
 * If the dstrank is DGRAM_BROADCAST, the transmission is a broadcast,
 * and should be delivered to all processors in the node.  If the dstrank
 * is DGRAM_ACKNOWLEDGE, the datagram is an acknowledgement datagram, in
 * which case the srcpe is the number of the acknowledger, the seqno is
 * always zero, and the user data is a list of the seqno's being
 * acknowledged.  There may be other dstrank codes for special functions.
 *
 * To send a message, one chops it up into datagrams and stores those
 * datagrams in a send-queue.  These outgoing datagrams aren't stored
 * in the explicit format shown above.  Instead, they are stored as
 * ImplicitDgrams, which contain the datagram header and a pointer to
 * the user data (which is in the user message buffer, which is in the
 * OutgoingMsg).  At transmission time these are combined together.

 * The combination of the datagram header with the user's data is
 * performed right in the user's message buffer.  Note that the
 * datagram header is exactly 64 bits.  One simply overwrites 64 bits
 * of the user's message with a datagram header, sends the datagram
 * straight from the user's message buffer, then restores the user's
 * buffer to its original state.  There is a small problem with the
 * first datagram of the message: one needs 64 bits of space to store
 * the datagram header.  To make sure this space is there, we added a
 * 64-bit unused space to the front of the Cmi message header.  In
 * addition to this, we also add 32 bits to the Cmi message header
 * to make room for a length-field, making it possible to identify
 * message boundaries.
 *
 * CONCURRENCY CONTROL
 *
 * This has changed recently.
 *
 * EFFICIENCY NOTES
 *
 * The sender-side does little copying.  The async and freeing send
 * routines do no copying at all.  The sync send routines copy the
 * message, then use the freeing-send routines.  The other alternative
 * is to not copy the message, and use the async send mechanism
 * combined with a blocking wait.  Blocking wait seems like a bad
 * idea, since it could take a VERY long time to get all those
 * datagrams out the door.
 *
 * The receiver side, unfortunately, must copy.  To avoid copying,
 * it would have to receive directly into a preallocated message buffer.
 * Unfortunately, this can't work: there's no way to know how much
 * memory to preallocate, and there's no way to know which datagram
 * is coming next.  Thus, we receive into fixed-size (large) datagram
 * buffers.  These are then inspected, and the messages extracted from
 * them.
 *
 * Note that we are allocating a large number of structs: OutgoingMsg's,
 * ImplicitDgrams, ExplicitDgrams.  By design, each of these structs
 * is a fixed-size structure.  Thus, we can do memory allocation by
 * simply keeping a linked-list of unused structs around.  The only
 * place where expensive memory allocation is performed is in the
 * sync routines.
 *
 * Since the datagrams from one node to another are fully ordered,
 * there is slightly more ordering than is needed: in theory, the
 * datagrams of one message don't need to be ordered relative to the
 * datagrams of another.  This was done to simplify the sequencing
 * mechanisms: implementing a fully-ordered stream is much simpler
 * than a partially-ordered one.  It also makes it possible to
 * modularize, layering the message transmitter on top of the
 * datagram-sequencer.  In other words, it was just easier this way.
 * Hopefully, this won't cause serious degradation: LAN's rarely get
 * datagrams out of order anyway.
 *
 * A potential efficiency problem is the lack of message-combining.
 * One datagram could conceivably contain several messages.  This
 * might be more efficient, it's not clear how much overhead is
 * involved in sending a short datagram.  Message-combining isn't
 * really ``integrated'' into the design of this software, but you
 * could fudge it as follows.  Whenever you pull a short datagram from
 * the send-queue, check the next one to see if it's also a short
 * datagram.  If so, pack them together into a ``combined'' datagram.
 * At the receive side, simply check for ``combined'' datagrams, and
 * treat them as if they were simply two datagrams.  This would
 * require extra copying.  I have no idea if this would be worthwhile.
 *
 *****************************************************************************/

/**
 * @addtogroup NET
 * @{
 */

/*****************************************************************************
 *
 * Include Files
 *
 ****************************************************************************/

#include <stdarg.h> /*<- was <varargs.h>*/

#define CMK_USE_PRINTF_HACK 0
#if CMK_USE_PRINTF_HACK
/*HACK: turn printf into CmiPrintf, by just defining our own
external symbol "printf".  This may be more trouble than it's worth,
since the only advantage is that it works properly with +syncprint.

This version *won't* work with fprintf(stdout,...) or C++ or Fortran I/O,
because they don't call printf.  Has to be defined up here because we probably 
haven't properly guessed this compiler's prototype for "printf".
*/
static void InternalPrintf(const char *f, va_list l);
int printf(const char *fmt, ...) {
	int nChar;
	va_list p; va_start(p, fmt);
        InternalPrintf(fmt,p);
	va_end(p);
	return 10;
}
#endif


#include "converse.h"
#include "memory-isomalloc.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>
#include <errno.h>
#include <setjmp.h>
#include <signal.h>
#include <string.h>

/* define machine debug */
#include "machine.h"

/******************* Producer-Consumer Queues ************************/
#include "pcqueue.h"

#include "machine-smp.h"

#if CMK_USE_POLL
#include <poll.h>
#endif

#if CMK_USE_GM
#include "gm.h"
struct gm_port *gmport = NULL;
int  portFinish = 0;
#endif

#include "conv-ccs.h"
#include "ccs-server.h"
#include "sockRoutines.h"

#if defined(_WIN32) && ! defined(__CYGWIN__)
/*For windows systems:*/
#  include <windows.h>
#  include <wincon.h>
#  include <sys/types.h>
#  include <sys/timeb.h>
#  define fdopen _fdopen
#  define SIGBUS -1  /*These signals don't exist in Win32*/
#  define SIGKILL -1
#  define SIGQUIT -1
#  define SIGTERM -1

#else /*UNIX*/
#  include <pwd.h>
#  include <unistd.h>
#  include <fcntl.h>
#  include <sys/file.h>
#endif

#if CMK_PERSISTENT_COMM
#include "persist_impl.h"
#endif

#define PRINTBUFSIZE 16384

/*
0: from smp thread
1: from interrupt
2: from worker thread
*/
static void CommunicationServer(int withDelayMs, int where);

void CmiHandleImmediate();
extern int CmemInsideMem();
extern void CmemCallWhenMemAvail();
static void ConverseRunPE(int everReturn);
void CmiYield(void);
void ConverseCommonExit(void);

static unsigned int dataport=0;
static int Cmi_mach_id=0; /* Machine-specific identifier (GM-only) */
static SOCKET       dataskt;

#define BROADCAST_SPANNING_FACTOR		4

/****************************************************************************
 *
 * Handling Errors
 *
 * Errors should be handled by printing a message on stderr and
 * calling exit(1).  Nothing should be sent to charmrun, no attempt at
 * communication should be made.  The other processes will notice the
 * abnormal termination and will deal with it.
 *
 * Rationale: if an error triggers an attempt to send a message,
 * the attempt to send a message is likely to trigger another error,
 * leading to an infinite loop and a process that spins instead of
 * shutting down.
 *
 *****************************************************************************/

static int machine_initiated_shutdown=0;
static int already_in_signal_handler=0;

static void CmiDestoryLocks();

static void machine_exit(int status)
{
  machine_initiated_shutdown=1;

  CmiDestoryLocks();		/* destory locks to prevent dead locking */

#if CMK_USE_GM
  if (gmport) { 
    gm_close(gmport); gmport = 0;
    gm_finalize();
  }
#endif
  exit(status);
}

static void charmrun_abort(const char*);

static void KillEveryone(const char *msg)
{
  charmrun_abort(msg);
  machine_exit(1);
}

static void KillEveryoneCode(n)
int n;
{
  char _s[100];
  sprintf(_s, "[%d] Fatal error #%d\n", CmiMyPe(), n);
  charmrun_abort(_s);
  machine_exit(1);
}

static void KillOnAllSigs(int sigNo)
{
  const char *sig="unknown signal";
  const char *suggestion="";
  if (machine_initiated_shutdown ||
      already_in_signal_handler) 
  	machine_exit(1); /*Don't infinite loop if there's a signal during a signal handler-- just die.*/
  already_in_signal_handler=1;

  CmiDestoryLocks();		/* destory locks */

  if (sigNo==SIGSEGV) {
     sig="segmentation violation";
     suggestion="Try running with '++debug', or linking with '-memory paranoid'.\n";
  }
  if (sigNo==SIGFPE) {
     sig="floating point exception";
     suggestion="Check for integer or floating-point division by zero.\n";
  }
  if (sigNo==SIGBUS) {
     sig="bus error";
     suggestion="Check for misaligned reads or writes to memory.\n";
  }
  if (sigNo==SIGILL) {
     sig="illegal instruction";
     suggestion="Check for calls to uninitialized function pointers.\n";
  }
  if (sigNo==SIGKILL) sig="caught signal KILL";
  if (sigNo==SIGQUIT) sig="caught signal QUIT";
  if (sigNo==SIGTERM) sig="caught signal TERM";
  MACHSTATE1(5,"     Caught signal %s ",sig);
/*ifdef this part*/
#ifdef __FAULT__
  if(sigNo == SIGKILL || sigNo == SIGQUIT || sigNo == SIGTERM){
		CmiPrintf("[%d] Caught but ignoring signal\n",CmiMyPe());
  }else{
#else
	{
#endif
   CmiError("------------- Processor %d Exiting: Caught Signal ------------\n"
  	"Signal: %s\n",CmiMyPe(),sig);
  	if (0!=suggestion[0])
  		CmiError("Suggestion: %s",suggestion);
  	CmiPrintStackTrace(1);
  	charmrun_abort(sig);
  	machine_exit(1);		
	}	
}

static void machine_atexit_check(void)
{
  if (!machine_initiated_shutdown)
    CmiAbort("unexpected call to exit by user program. Must use CkExit, not exit!");
  printf("Program finished.\n");
#if 0 /*Wait for the user to press any key (for Win32 debugging)*/
  fgetc(stdin);
#endif
}

#if !defined(_WIN32) || defined(__CYGWIN__)
static void HandleUserSignals(int signum)
{
  int condnum = ((signum==SIGUSR1) ? CcdSIGUSR1 : CcdSIGUSR2);
  CcdRaiseCondition(condnum);
}
#endif

static void PerrorExit(const char *msg)
{
  perror(msg);
  machine_exit(1);
}

/*****************************************************************************
 *
 *     Utility routines for network machine interface.
 *
 *****************************************************************************/

/*
Horrific #defines to hide the differences between select() and poll().
 */
#if CMK_USE_POLL /*poll() version*/
# define CMK_PIPE_DECL(delayMs) \
	struct pollfd fds[10]; \
	int nFds_sto=0; int *nFds=&nFds_sto; \
	int pollDelayMs=delayMs;
# define CMK_PIPE_SUB fds,nFds
# define CMK_PIPE_CALL() poll(fds, *nFds, pollDelayMs); *nFds=0

# define CMK_PIPE_PARAM struct pollfd *fds,int *nFds
# define CMK_PIPE_ADDREAD(rd_fd) \
	do {fds[*nFds].fd=rd_fd; fds[*nFds].events=POLLIN; (*nFds)++;} while(0)
# define CMK_PIPE_ADDWRITE(wr_fd) \
	do {fds[*nFds].fd=wr_fd; fds[*nFds].events=POLLOUT; (*nFds)++;} while(0)
# define CMK_PIPE_CHECKREAD(rd_fd) fds[(*nFds)++].revents&POLLIN
# define CMK_PIPE_CHECKWRITE(wr_fd) fds[(*nFds)++].revents&POLLOUT

#else /*select() version*/

# define CMK_PIPE_DECL(delayMs) \
	fd_set rfds_sto,wfds_sto;\
	fd_set *rfds=&rfds_sto,*wfds=&wfds_sto; struct timeval tmo; \
	FD_ZERO(rfds); FD_ZERO(wfds);tmo.tv_sec=0; tmo.tv_usec=1000*delayMs;
# define CMK_PIPE_SUB rfds,wfds
# define CMK_PIPE_CALL() select(FD_SETSIZE, rfds, wfds, NULL, &tmo)

# define CMK_PIPE_PARAM fd_set *rfds,fd_set *wfds
# define CMK_PIPE_ADDREAD(rd_fd) FD_SET(rd_fd,rfds)
# define CMK_PIPE_ADDWRITE(wr_fd) FD_SET(wr_fd,wfds)
# define CMK_PIPE_CHECKREAD(rd_fd) FD_ISSET(rd_fd,rfds)
# define CMK_PIPE_CHECKWRITE(wr_fd) FD_ISSET(wr_fd,wfds)
#endif

static void CMK_PIPE_CHECKERR(void) {
#if defined(_WIN32) && !defined(__CYGWIN__)
/* Win32 socket seems to randomly return inexplicable errors
here-- WSAEINVAL, WSAENOTSOCK-- yet everything is actually OK. 
        int err=WSAGetLastError();
        CmiPrintf("(%d)Select returns -1; errno=%d, WSAerr=%d\n",withDelayMs,errn
o,err);
*/
#else /*UNIX machine*/
        if (errno!=EINTR)
                KillEveryone("Socket error in CheckSocketsReady!\n");
#endif
}


static void CmiStdoutFlush(void);
static int  CmiStdoutNeedsService(void);
static void CmiStdoutService(void);
static void CmiStdoutAdd(CMK_PIPE_PARAM);
static void CmiStdoutCheck(CMK_PIPE_PARAM);


double GetClock(void)
{
#if defined(_WIN32) && !defined(__CYGWIN__)
  struct _timeb tv; 
  _ftime(&tv);
  return (tv.time * 1.0 + tv.millitm * 1.0E-3);
#else
  struct timeval tv; int ok;
  ok = gettimeofday(&tv, NULL);
  if (ok<0) { perror("gettimeofday"); KillEveryoneCode(9343112); }
  return (tv.tv_sec * 1.0 + tv.tv_usec * 1.0E-6);
#endif
}

char *CopyMsg(char *msg, int len)
{
  char *copy = (char *)CmiAlloc(len);
  if (!copy)
      fprintf(stderr, "Out of memory\n");
  memcpy(copy, msg, len);
  return copy;
}

/***********************************************************************
 *
 * Abort function:
 *
 ************************************************************************/

static int  Cmi_truecrash;
static int already_aborting=0;
void CmiAbort(const char *message)
{
  if (already_aborting) machine_exit(1);
  already_aborting=1;
  MACHSTATE1(5,"CmiAbort(%s)",message);
  
  CmiError("------------- Processor %d Exiting: Called CmiAbort ------------\n"
  	"Reason: %s\n",CmiMyPe(),message);
  CmiPrintStackTrace(0);
  
  /*Send off any remaining prints*/
  CmiStdoutFlush();
  
  if(Cmi_truecrash) {
    printf("CHARM++ FATAL ERROR: %s\n", message);
    *(int *)NULL = 0; /*Write to null, causing bus error*/
  } else {
    charmrun_abort(message);
    machine_exit(1);
  }
}


/******************************************************************************
 *
 * CmiEnableAsyncIO
 *
 * The net and tcp versions use a bunch of unix processes talking to each
 * other via file descriptors.  We need for a signal SIGIO to be generated
 * each time a message arrives, making it possible to write a signal
 * handler to handle the messages.  The vast majority of unixes can,
 * in fact, do this.  However, there isn't any standard for how this is
 * supposed to be done, so each version of UNIX has a different set of
 * calls to turn this signal on.  So, there is like one version here for
 * every major brand of UNIX.
 *
 *****************************************************************************/

#if CMK_ASYNC_USE_F_SETFL_AND_F_SETOWN
#include <fcntl.h>
void CmiEnableAsyncIO(int fd)
{
  if ( fcntl(fd, F_SETOWN, getpid()) < 0 ) {
    CmiError("setting socket owner: %s\n", strerror(errno)) ;
    exit(1);
  }
  if ( fcntl(fd, F_SETFL, FASYNC) < 0 ) {
    CmiError("setting socket async: %s\n", strerror(errno)) ;
    exit(1);
  }
}
#else
void CmiEnableAsyncIO(int fd) { }
#endif

/* We should probably have a set of "CMK_NONBLOCK_USE_..." defines here:*/
#if !defined(_WIN32) || defined(__CYGWIN__)
void CmiEnableNonblockingIO(int fd) {
  int on=1;
  if (fcntl(fd,F_SETFL,O_NONBLOCK,&on)<0) {
    CmiError("setting nonblocking IO: %s\n", strerror(errno)) ;
    exit(1);
  }
}
#else
void CmiEnableNonblockingIO(int fd) { }
#endif


/******************************************************************************
 *
 * Configuration Data
 *
 * This data is all read in from the NETSTART variable (provided by the
 * charmrun) and from the command-line arguments.  Once read in, it is never
 * modified.
 *
 *****************************************************************************/

int               _Cmi_mynode;    /* Which address space am I */
int               _Cmi_mynodesize;/* Number of processors in my address space */
int               _Cmi_numnodes;  /* Total number of address spaces */
int               _Cmi_numpes;    /* Total number of processors */
static int        Cmi_nodestart; /* First processor in this address space */
static skt_ip_t   Cmi_self_IP;
static skt_ip_t   Cmi_charmrun_IP; /*Address of charmrun machine*/
static int        Cmi_charmrun_port;
static int        Cmi_charmrun_pid;
static int        Cmi_charmrun_fd=-1;

static int    Cmi_netpoll;
static int    Cmi_asyncio;
static int    Cmi_idlepoll;
static int    Cmi_syncprint;
static int Cmi_print_stats = 0;

static void parse_netstart(void)
{
  char *ns;
  int nread;
  int port;
  ns = getenv("NETSTART");
  if (ns!=0) 
  {/*Read values set by Charmrun*/
        char Cmi_charmrun_name[1024];
        nread = sscanf(ns, "%d%s%d%d%d",
                 &_Cmi_mynode,
                 Cmi_charmrun_name, &Cmi_charmrun_port,
                 &Cmi_charmrun_pid, &port);
	Cmi_charmrun_IP=skt_lookup_ip(Cmi_charmrun_name);

        if (nread!=5) {
                fprintf(stderr,"Error parsing NETSTART '%s'\n",ns);
                exit(1);
        }
  } else 
  {/*No charmrun-- set flag values for standalone operation*/
  	_Cmi_mynode=0;
  	Cmi_charmrun_IP=_skt_invalid_ip;
  	Cmi_charmrun_port=0;
  	Cmi_charmrun_pid=0;
        dataport = -1;
  }
}

static void extract_common_args(char **argv)
{
  if (CmiGetArgFlagDesc(argv,"+stats","Print network statistics at shutdown"))
    Cmi_print_stats = 1;
}

/* for SMP */
#include "machine-smp.c"

CsvDeclare(CmiNodeState, NodeState);

/* Immediate message support */
#define CMI_DEST_RANK(msg)	*(int *)(msg)
#include "immediate.c"

/******************************************************************************
 *
 * Packet Performance Logging
 *
 * This module is designed to give a detailed log of the packets and their
 * acknowledgements, for performance tuning.  It can be disabled.
 *
 *****************************************************************************/

#define LOGGING 0

#if LOGGING

typedef struct logent {
  double time;
  int seqno;
  int srcpe;
  int dstpe;
  int kind;
} *logent;


logent log;
int    log_pos;
int    log_wrap;

static void log_init(void)
{
  log = (logent)malloc(50000 * sizeof(struct logent));
  _MEMCHECK(log);
  log_pos = 0;
  log_wrap = 0;
}

static void log_done(void)
{
  char logname[100]; FILE *f; int i, size;
  sprintf(logname, "log.%d", _Cmi_mynode);
  f = fopen(logname, "w");
  if (f==0) KillEveryone("fopen problem");
  if (log_wrap) size = 50000; else size=log_pos;
  for (i=0; i<size; i++) {
    logent ent = log+i;
    fprintf(f, "%1.4f %d %c %d %d\n",
	    ent->time, ent->srcpe, ent->kind, ent->dstpe, ent->seqno);
  }
  fclose(f);
}

void printLog(void)
{
  char logname[100]; FILE *f; int i, j, size;
  static int logged = 0;
  if (logged)
      return;
  logged = 1;
  CmiPrintf("Logging: %d\n", _Cmi_mynode);
  sprintf(logname, "log.%d", _Cmi_mynode);
  f = fopen(logname, "w");
  if (f==0) KillEveryone("fopen problem");
  for (i = 5000; i; i--)
  {
  /*for (i=0; i<size; i++) */
    j = log_pos - i;
    if (j < 0)
    {
        if (log_wrap)
	    j = 5000 + j;
	else
	    j = 0;
    };
    {
    logent ent = log+j;
    fprintf(f, "%1.4f %d %c %d %d\n",
	    ent->time, ent->srcpe, ent->kind, ent->dstpe, ent->seqno);
    }
  }
  fclose(f);
  CmiPrintf("Done Logging: %d\n", _Cmi_mynode);
}

#define LOG(t,s,k,d,q) { if (log_pos==50000) { log_pos=0; log_wrap=1;} { logent ent=log+log_pos; ent->time=t; ent->srcpe=s; ent->kind=k; ent->dstpe=d; ent->seqno=q; log_pos++; }}

#endif

#if !LOGGING

#define log_init() /*empty*/
#define log_done() /*empty*/
#define printLog() /*empty*/
#define LOG(t,s,k,d,q) /*empty*/

#endif

/******************************************************************************
 *
 * Node state
 *
 *****************************************************************************/


static CmiNodeLock    Cmi_scanf_mutex;
static double         Cmi_clock;
static double         Cmi_check_delay = 3.0;

/******************************************************************************
 *
 * OS Threads
 * SMP implementation moved to machine-smp.c
 *****************************************************************************/

/************************ No kernel SMP threads ***************/
#if CMK_SHARED_VARS_UNAVAILABLE

static volatile int memflag=0;
void CmiMemLock() { memflag++; }
void CmiMemUnlock() { memflag--; }

static volatile int comm_flag=0;
#define CmiCommLockOrElse(dothis) if (comm_flag!=0) dothis
#ifndef MACHLOCK_DEBUG
#  define CmiCommLock() (comm_flag=1)
#  define CmiCommUnlock() (comm_flag=0)
#else /* Error-checking flag locks */
void CmiCommLock(void) {
  MACHLOCK_ASSERT(!comm_flag,"CmiCommLock");
  comm_flag=1;
}
void CmiCommUnlock(void) {
  MACHLOCK_ASSERT(comm_flag,"CmiCommUnlock");
  comm_flag=0;
}
#endif

static struct CmiStateStruct Cmi_state;
int _Cmi_mype;
int _Cmi_myrank=0; /* Normally zero; only 1 during SIGIO handling */
#define CmiGetState() (&Cmi_state)
#define CmiGetStateN(n) (&Cmi_state)

void CmiYield(void) { sleep(0); }

static void CommunicationInterrupt(int ignored)
{
  MACHLOCK_ASSERT(!_Cmi_myrank,"CommunicationInterrupt");
  if (memflag || comm_flag || _immRunning || CmiCheckImmediateLock(0)) 
  { /* Already busy inside malloc, comm, or immediate messages */
    MACHSTATE(5,"--SKIPPING SIGIO--");
    return;
  }
  MACHSTATE1(2,"--BEGIN SIGIO comm_flag: %d--", comm_flag)
  {
    /*Make sure any malloc's we do in here are NOT migratable:*/
    CmiIsomallocBlockList *oldList=CmiIsomallocBlockListActivate(NULL);
/*    _Cmi_myrank=1; */
    CommunicationServer(0, 1);     /* from interrupt */
/*    _Cmi_myrank=0; */
    CmiIsomallocBlockListActivate(oldList);
  }
  MACHSTATE(2,"--END SIGIO--")
}

extern void CmiSignal(int sig1, int sig2, int sig3, void (*handler)());

static void CmiStartThreads(char **argv)
{
  if ((_Cmi_numpes != _Cmi_numnodes) || (_Cmi_mynodesize != 1))
    KillEveryone
      ("Multiple cpus unavailable, don't use cpus directive in nodesfile.\n");
  
  CmiStateInit(Cmi_nodestart, 0, &Cmi_state);
  _Cmi_mype = Cmi_nodestart;

  /* Prepare Cpv's for immediate messages: */
  _Cmi_myrank=1;
  CommunicationServerInit();
  _Cmi_myrank=0;
  
#if !CMK_ASYNC_NOT_NEEDED
  if (Cmi_asyncio)
  {
    CmiSignal(SIGIO, 0, 0, CommunicationInterrupt);
    if (!Cmi_netpoll) {
      if (dataskt!=-1) CmiEnableAsyncIO(dataskt);
      if (Cmi_charmrun_fd!=-1) CmiEnableAsyncIO(Cmi_charmrun_fd);
    }
  }
#endif
}

static void CmiDestoryLocks()
{
  comm_flag = 0;
  memflag = 0;
}

#endif

/*Network progress utility variables. Period controls the rate at
  which the network poll is called */

CpvDeclare(int , networkProgressCount);
int networkProgressPeriod;

CpvDeclare(void *, CmiLocalQueue);


#ifndef CmiMyPe
int CmiMyPe(void) 
{ 
  return CmiGetState()->pe; 
}
#endif
#ifndef CmiMyRank
int CmiMyRank(void)
{
  return CmiGetState()->rank;
}
#endif

/*Add a message to this processor's receive queue 
  Must be called while holding comm. lock
*/
static void CmiPushPE(int pe,void *msg)
{
  CmiState cs=CmiGetStateN(pe);
  MACHSTATE1(2,"Pushing message into %d's queue",pe);  
  MACHLOCK_ASSERT(comm_flag,"CmiPushPE")

#if CMK_IMMEDIATE_MSG
  if (CmiIsImmediate(msg)) {
    CmiPushImmediateMsg(msg);
    return;
  }
#endif
  PCQueuePush(cs->recv,msg);
  CmiIdleLock_addMessage(&cs->idle);
}

#if CMK_NODE_QUEUE_AVAILABLE
/*Add a message to the node queue.  
  Must be called while holding comm. lock
*/
static void CmiPushNode(void *msg)
{
  CmiState cs=CmiGetStateN(0);
  MACHSTATE(2,"Pushing message into node queue");
  MACHLOCK_ASSERT(comm_flag,"CmiPushNode")
  
#if CMK_IMMEDIATE_MSG
  if (CmiIsImmediate(msg)) {
    MACHSTATE(2,"Pushing Immediate message into queue");
    CmiPushImmediateMsg(msg);
    return;
  }
#endif 
  PCQueuePush(CsvAccess(NodeState).NodeRecv,msg);
  /*Silly: always try to wake up processor 0, so at least *somebody*
    will be awake to handle the message*/
  CmiIdleLock_addMessage(&cs->idle);
}
#endif

/***************************************************************
 Communication with charmrun:
 We can send (ctrl_sendone) and receive (ctrl_getone)
 messages on a TCP socket connected to charmrun.
 This is used for printfs, CCS, etc; and also for
 killing ourselves if charmrun dies.
*/

/*This flag prevents simultanious outgoing
messages on the charmrun socket.  It is protected
by the commlock.*/
static int Cmi_charmrun_fd_sendflag=0;

/* ctrl_sendone */
static int sendone_abort_fn(int code,const char *msg) {
	fprintf(stderr,"Socket error %d in ctrl_sendone! %s\n",code,msg);
	machine_exit(1);
	return -1;
}

static void ctrl_sendone_nolock(const char *type,
				const char *data1,int dataLen1,
				const char *data2,int dataLen2)
{
  const void *bufs[3]; int lens[3]; int nBuffers=0;
  ChMessageHeader hdr;
  skt_abortFn oldAbort=skt_set_abort(sendone_abort_fn);
  if (Cmi_charmrun_fd==-1) 
  	charmrun_abort("ctrl_sendone called in standalone!\n");
  Cmi_charmrun_fd_sendflag=1;
  ChMessageHeader_new(type,dataLen1+dataLen2,&hdr);
  bufs[nBuffers]=&hdr; lens[nBuffers]=sizeof(hdr); nBuffers++;
  if (dataLen1>0) {bufs[nBuffers]=data1; lens[nBuffers]=dataLen1; nBuffers++;}
  if (dataLen2>0) {bufs[nBuffers]=data2; lens[nBuffers]=dataLen2; nBuffers++;}
  skt_sendV(Cmi_charmrun_fd,nBuffers,bufs,lens);
  Cmi_charmrun_fd_sendflag=0;
  skt_set_abort(oldAbort);
}

static void ctrl_sendone_locking(const char *type,
				const char *data1,int dataLen1,
				const char *data2,int dataLen2)
{
  CmiCommLock();
  ctrl_sendone_nolock(type,data1,dataLen1,data2,dataLen2);
  CmiCommUnlock();
}

static double Cmi_check_last;

/* if charmrun dies, we finish */
static void pingCharmrun(void *ignored) 
{
  double clock=GetClock();
  if (clock > Cmi_check_last + Cmi_check_delay) {
    MACHSTATE2(2,"CommunicationsClock pinging charmrun comm_flag=%d Cmi_charmrun_fd_sendflag=%d", comm_flag, Cmi_charmrun_fd_sendflag);
    Cmi_check_last = clock; 
#if CMK_USE_GM
    if (!Cmi_netpoll)  /* GM netpoll, charmrun service is done in interrupt */
#endif
    CmiCommLockOrElse(return;); /*Already busy doing communication*/
    if (Cmi_charmrun_fd_sendflag) return; /*Busy talking to charmrun*/
    CmiCommLock();
    ctrl_sendone_nolock("ping",NULL,0,NULL,0); /*Charmrun may have died*/
    CmiCommUnlock();
  }
#if 1
#if CMK_USE_GM
  if (!Cmi_netpoll)
#endif
  CmiStdoutFlush(); /*Make sure stdout buffer hasn't filled up*/
#endif
}

/* periodic charm ping, for gm and netpoll */
static void pingCharmrunPeriodic(void *ignored)
{
  pingCharmrun(ignored);
  CcdCallFnAfter((CcdVoidFn)pingCharmrunPeriodic,NULL,1000);
}

static int ignore_further_errors(int c,const char *msg) {machine_exit(2);return -1;}
static void charmrun_abort(const char *s)
{
  if (Cmi_charmrun_fd==-1) {/*Standalone*/
  	fprintf(stderr,"Charm++ fatal error:\n%s\n",s);
  	abort();
  } else {
	char msgBuf[80];
  	skt_set_abort(ignore_further_errors);
	sprintf(msgBuf,"Fatal error on PE %d> ",CmiMyPe());
  	ctrl_sendone_nolock("abort",msgBuf,strlen(msgBuf),s,strlen(s)+1);
  }
}

/* ctrl_getone */

#ifdef __FAULT__
#include "machine-recover.c"
#endif

static void node_addresses_store(ChMessage *msg);

static void ctrl_getone(void)
{
  ChMessage msg;
  MACHSTATE(2,"ctrl_getone")
  MACHLOCK_ASSERT(comm_flag,"ctrl_getone")
  ChMessage_recv(Cmi_charmrun_fd,&msg);
  MACHSTATE1(2,"ctrl_getone recv one '%s'", msg.header.type)

  if (strcmp(msg.header.type,"die")==0) {
    MACHSTATE(2,"ctrl_getone bye bye")
    fprintf(stderr,"aborting: %s\n",msg.data);
    log_done();
    ConverseCommonExit();
    machine_exit(0);
  } 
#if CMK_CCS_AVAILABLE
  else if (strcmp(msg.header.type, "req_fw")==0) {
    CcsImplHeader *hdr=(CcsImplHeader *)msg.data;
	/*Sadly, I *can't* do a:
      CcsImpl_netRequest(hdr,msg.data+sizeof(CcsImplHeader));
	here, because I can't send converse messages in the
	communication thread.  I *can* poke this message into 
	any convenient processor's queue, though:  (OSL, 9/14/2000)
	*/
    int pe=0;/*<- node-local processor number. Any one will do.*/
    void *cmsg=(void *)CcsImpl_ccs2converse(hdr,msg.data+sizeof(CcsImplHeader),NULL);
    MACHSTATE(2,"Incoming CCS request");
    CmiPushPE(pe,cmsg);
  }
#endif
#ifdef __FAULT__	
  else if(strcmp(msg.header.type,"crashnode")==0) {
	crash_node_handle(&msg);
  }
  else if(strcmp(msg.header.type,"initnodetab")==0) {
	/** A processor crashed and got recreated. So charmrun sent 
	  across the whole nodetable data to update this processor*/
	node_addresses_store(&msg);
	fprintf(stdout,"nodetable added %d\n",CmiMyPe());
  }
#endif
  else {
  /* We do not use KillEveryOne here because it calls CmiMyPe(),
   * which is not available to the communication thread on an SMP version.
   */
    charmrun_abort("ERROR> Unrecognized message from charmrun.\n");
    machine_exit(1);
  }
  
  MACHSTATE(2,"ctrl_getone done")
  ChMessage_free(&msg);
}

#if CMK_CCS_AVAILABLE
/*Deliver this reply data to this reply socket.
  The data is forwarded to CCS server via charmrun.*/
void CcsImpl_reply(CcsImplHeader *hdr,int repLen,const void *repData)
{
  MACHSTATE(2,"Outgoing CCS reply");
  ctrl_sendone_locking("reply_fw",(const char *)hdr,sizeof(CcsImplHeader),
		       repData,repLen);
  MACHSTATE(1,"Outgoing CCS reply away");
}
#endif

/*****************************************************************************
 *
 * CmiPrintf, CmiError, CmiScanf
 *
 *****************************************************************************/
static void InternalWriteToTerminal(int isStdErr,const char *str,int len);
static void InternalPrintf(const char *f, va_list l)
{
  ChMessage replymsg;
  char *buffer = CmiTmpAlloc(PRINTBUFSIZE);
  CmiStdoutFlush();
  vsprintf(buffer, f, l);
  if(Cmi_syncprint) {
	  CmiCommLock();
  	  ctrl_sendone_nolock("printsyn", buffer,strlen(buffer)+1,NULL,0);
  	  ChMessage_recv(Cmi_charmrun_fd,&replymsg);
  	  ChMessage_free(&replymsg);
	  CmiCommUnlock();
  } else {
  	  ctrl_sendone_locking("print", buffer,strlen(buffer)+1,NULL,0);
  }
  InternalWriteToTerminal(0,buffer,strlen(buffer));
  CmiTmpFree(buffer);
}

static void InternalError(const char *f, va_list l)
{
  ChMessage replymsg;
  char *buffer = CmiTmpAlloc(PRINTBUFSIZE);
  CmiStdoutFlush();
  vsprintf(buffer, f, l);
  if(Cmi_syncprint) {
  	  ctrl_sendone_locking("printerrsyn", buffer,strlen(buffer)+1,NULL,0);
	  CmiCommLock();
  	  ChMessage_recv(Cmi_charmrun_fd,&replymsg);
  	  ChMessage_free(&replymsg);
	  CmiCommUnlock();
  } else {
  	  ctrl_sendone_locking("printerr", buffer,strlen(buffer)+1,NULL,0);
  }
  InternalWriteToTerminal(1,buffer,strlen(buffer));
  CmiTmpFree(buffer);
}

static int InternalScanf(char *fmt, va_list l)
{
  ChMessage replymsg;
  char *ptr[20];
  char *p; int nargs, i;
  nargs=0;
  p=fmt;
  while (*p) {
    if ((p[0]=='%')&&(p[1]=='*')) { p+=2; continue; }
    if ((p[0]=='%')&&(p[1]=='%')) { p+=2; continue; }
    if (p[0]=='%') { nargs++; p++; continue; }
    if (*p=='\n') *p=' '; p++;
  }
  if (nargs > 18) KillEveryone("CmiScanf only does 18 args.\n");
  for (i=0; i<nargs; i++) ptr[i]=va_arg(l, char *);
  CmiLock(Cmi_scanf_mutex);
  if (Cmi_charmrun_fd!=-1)
  {/*Send charmrun the format string*/
        ctrl_sendone_locking("scanf", fmt, strlen(fmt)+1,NULL,0);
        /*Wait for the reply (characters to scan) from charmrun*/
        CmiCommLock();
        ChMessage_recv(Cmi_charmrun_fd,&replymsg);
        i = sscanf((char*)replymsg.data, fmt,
                     ptr[ 0], ptr[ 1], ptr[ 2], ptr[ 3], ptr[ 4], ptr[ 5],
                     ptr[ 6], ptr[ 7], ptr[ 8], ptr[ 9], ptr[10], ptr[11],
                     ptr[12], ptr[13], ptr[14], ptr[15], ptr[16], ptr[17]);
        ChMessage_free(&replymsg);
        CmiCommUnlock();
  } else
  {/*Just do the scanf normally*/
        i=scanf(fmt, ptr[ 0], ptr[ 1], ptr[ 2], ptr[ 3], ptr[ 4], ptr[ 5],
                     ptr[ 6], ptr[ 7], ptr[ 8], ptr[ 9], ptr[10], ptr[11],
                     ptr[12], ptr[13], ptr[14], ptr[15], ptr[16], ptr[17]);
  }
  CmiUnlock(Cmi_scanf_mutex);
  return i;
}

#if CMK_CMIPRINTF_IS_A_BUILTIN

/*New stdarg.h declarations*/
void CmiPrintf(const char *fmt, ...)
{
  va_list p; va_start(p, fmt);
  if (Cmi_charmrun_fd!=-1)
    InternalPrintf(fmt, p);
  else
    vfprintf(stdout,fmt,p);
  va_end(p);
}

void CmiError(const char *fmt, ...)
{
  va_list p; va_start (p, fmt);
  if (Cmi_charmrun_fd!=-1)
    InternalError(fmt, p);
  else
    vfprintf(stderr,fmt,p);
  va_end(p);
}

int CmiScanf(const char *fmt, ...)
{
  va_list p; int i; va_start(p, fmt);
  i = InternalScanf((char *)fmt, p);
  va_end(p);
  return i;
}

#endif

/***************************************************************************
 * Output redirection:
 *  When people don't use CkPrintf, like above, we'd still like to be able
 * to collect their output.  Thus we make a pipe and dup2 it to stdout,
 * which lets us read the characters sent to stdout at our lesiure.
 ***************************************************************************/

/*Can read from stdout or stderr using these fd's*/
static int readStdout[2]; 
static int writeStdout[2]; /*The original stdout/stderr sockets*/ 
static int serviceStdout[2]; /*(bool) Normally zero; one if service needed.*/
#define readStdoutBufLen (16*1024)
static char readStdoutBuf[readStdoutBufLen+1]; /*Protected by comm. lock*/
static int servicingStdout;

/*Initialization-- should only be called once per node*/
static void CmiStdoutInit(void) {
	int i;
	if (Cmi_charmrun_fd==-1) return; /* standalone mode */
	
/*There's some way to do this same thing in windows, but I don't know how*/
#if !defined(_WIN32) || defined(__CYGWIN__)
	/*Prevent buffering in stdio library:*/
	setbuf(stdout,NULL); setbuf(stderr,NULL);

	/*Reopen stdout and stderr fd's as new pipes:*/
        for (i=0;i<2;i++) {
		int pair[2];
		int srcFd=1+i; /* 1 is stdout; 2 is stderr */
		
		/*First, save a copy of the original stdout*/
		writeStdout[i]=dup(srcFd);
#if 0
		/*Build a pipe to connect to stdout (4kb buffer, but no SIGIO...)*/
		if (-1==pipe(pair)) {perror("building stdio redirection pipe"); exit(1);}
#else
	       /* UNIX socket (16kb default buffer, and works with SIGIO!) */
		if (-1==socketpair(PF_UNIX,SOCK_STREAM,0,pair)) 
			{perror("building stdio redirection socketpair"); exit(1);}
#endif
		readStdout[i]=pair[0]; /*We get the read end of pipe*/
		if (-1==dup2(pair[1],srcFd)) {perror("dup2 redirection pipe"); exit(1);}
		
#if 0 /*Keep writes from blocking.  This just drops excess output, which is bad.*/
		CmiEnableNonblockingIO(srcFd);
#endif
#if CMK_SHARED_VARS_UNAVAILABLE
                if (Cmi_asyncio)
		{
  /*No communication thread-- get a SIGIO on each write(), which keeps the buffer clean*/
			CmiEnableAsyncIO(readStdout[i]);
		}
#endif
	}
#else
/*Windows system-- just fake reads for now*/
# ifndef read
#  define read(x,y,z) 0
# endif
# ifndef write
#  define write(x,y,z) 
# endif
#endif
}

/*Sends data to original stdout (e.g., for ++debug or ++in-xterm)*/
static void InternalWriteToTerminal(int isStdErr,const char *str,int len)
{
	write(writeStdout[isStdErr],str,len);	
}

/*
  Service this particular stdout pipe.  
  Must hold comm. lock.
*/
static void CmiStdoutServiceOne(int i) {
	int nBytes;
	const static char *cmdName[2]={"print","printerr"};
	servicingStdout=1;
	while(1) {
		const char *tooMuchWarn=NULL; int tooMuchLen=0;
		if (!skt_select1(readStdout[i],0)) break; /*Nothing to read*/
		nBytes=read(readStdout[i],readStdoutBuf,readStdoutBufLen);
		if (nBytes<=0) break; /*Nothing to send*/
		
		/*Send these bytes off to charmrun*/
		readStdoutBuf[nBytes]=0; /*Zero-terminate read string*/
		nBytes++; /*Include zero-terminator in message to charmrun*/
		
		if (nBytes>=readStdoutBufLen-100) 
		{ /*We must have filled up our output pipe-- most output libraries
		   don't handle this well (e.g., glibc printf just drops the line).*/
			
			tooMuchWarn="\nWARNING: Too much output at once-- possible output discontinuity!\n"
				"Use CkPrintf to avoid discontinuity (and this warning).\n\n";
			nBytes--; /*Remove terminator from user's data*/
			tooMuchLen=strlen(tooMuchWarn)+1;
		}
		ctrl_sendone_nolock(cmdName[i],readStdoutBuf,nBytes,
				    tooMuchWarn,tooMuchLen);
		
		InternalWriteToTerminal(i,readStdoutBuf,nBytes);
	}
	servicingStdout=0;
	serviceStdout[i]=0; /*This pipe is now serviced*/
}

/*Service all stdout pipes, whether it looks like they need it
  or not.  Used when you aren't sure if select() has been called recently.
  Must hold comm. lock.
*/
static void CmiStdoutServiceAll(void) {
	int i;
	for (i=0;i<2;i++) {
		if (readStdout[i]==0) continue; /*Pipe not open*/
		CmiStdoutServiceOne(i);
	}
}

/*Service any outstanding stdout pipes.
  Must hold comm. lock.
*/
static void CmiStdoutService(void) {
	CmiStdoutServiceAll();
}

/*Add our pipes to the pile for select() or poll().
  Both can be called with or without the comm. lock.
*/
static void CmiStdoutAdd(CMK_PIPE_PARAM) {
	int i;
	for (i=0;i<2;i++) {
		if (readStdout[i]==0) continue; /*Pipe not open*/
		CMK_PIPE_ADDREAD(readStdout[i]);
	}
}
static void CmiStdoutCheck(CMK_PIPE_PARAM) {
	int i;
	for (i=0;i<2;i++) {
		if (readStdout[i]==0) continue; /*Pipe not open*/
		if (CMK_PIPE_CHECKREAD(readStdout[i])) serviceStdout[i]=1;
	}
}
static int CmiStdoutNeedsService(void) {
	return (serviceStdout[0]!=0 || serviceStdout[1]!=0);
}

/*Called every few milliseconds to flush the stdout pipes*/
static void CmiStdoutFlush(void) {
	if (servicingStdout) return; /* might be called by SIGALRM */
	CmiCommLockOrElse( return; )
	CmiCommLock();
	CmiStdoutServiceAll();
	CmiCommUnlock();
}

/***************************************************************************
 * Message Delivery:
 *
 ***************************************************************************/

#include "machine-dgram.c"


#ifndef CmiNodeFirst
int CmiNodeFirst(int node) { return nodes[node].nodestart; }
int CmiNodeSize(int node)  { return nodes[node].nodesize; }
#endif

#ifndef CmiNodeOf
int CmiNodeOf(int pe)      { return (nodes_by_pe[pe] - nodes); }
int CmiRankOf(int pe)      { return pe - (nodes_by_pe[pe]->nodestart); }
#endif


/*****************************************************************************
 *
 * node_addresses
 *
 *  These two functions fill the node-table.
 *
 *
 *   This node, like all others, first sends its own address to charmrun
 *   using this command:
 *
 *     Type: nodeinfo
 *     Data: Big-endian 4-byte ints
 *           <my-node #><Dataport>
 *
 *   When charmrun has all the addresses, he sends this table to me:
 *
 *     Type: nodes
 *     Data: Big-endian 4-byte ints
 *           <number of nodes n>
 *           <#PEs><IP><Dataport> Node 0
 *           <#PEs><IP><Dataport> Node 1
 *           ...
 *           <#PEs><IP><Dataport> Node n-1
 *
 *****************************************************************************/

/*Note: node_addresses_obtain is called before starting
  threads, so no locks are needed (or valid!)*/
static void node_addresses_obtain(char **argv)
{
  ChMessage nodetabmsg; /* info about all nodes*/
  if (Cmi_charmrun_fd==-1) 
  {/*Standalone-- fake a single-node nodetab message*/
  	int npes=1;
  	ChSingleNodeinfo *fakeTab;
	ChMessage_new("nodeinfo",sizeof(ChSingleNodeinfo),&nodetabmsg);
	fakeTab=(ChSingleNodeinfo *)(nodetabmsg.data);
  	CmiGetArgIntDesc(argv,"+p",&npes,"Set the number of processes to create");
#if CMK_SHARED_VARS_UNAVAILABLE
	if (npes!=1) {
		fprintf(stderr,
			"To use multiple processors, you must run this program as:\n"
			" > charmrun +p%d %s <args>\n"
			"or build the %s-smp version of Charm++.\n",
			npes,argv[0],CMK_MACHINE_NAME);
		exit(1);
	}
#endif
	/*This is a stupid hack: we expect the *number* of nodes
	followed by ChNodeinfo structs; so we use a ChSingleNodeinfo
	(which happens to have exactly that layout!) and stuff
	a 1 into the "node number" slot
	*/
	fakeTab->nodeNo=ChMessageInt_new(1); /* <- hack */
	fakeTab->info.nPE=ChMessageInt_new(npes);
	fakeTab->info.dataport=ChMessageInt_new(0);
	fakeTab->info.IP=_skt_invalid_ip;
  }
  else 
  { /*Contact charmrun for machine info.*/
	ChSingleNodeinfo me;

  	me.nodeNo=ChMessageInt_new(_Cmi_mynode);
	/*The nPE and IP fields are set by charmrun--
	  these values don't matter.
	*/
	me.info.nPE=ChMessageInt_new(0);
	me.info.IP=_skt_invalid_ip;
	me.info.mach_id=ChMessageInt_new(Cmi_mach_id);
  	me.info.dataport=ChMessageInt_new(dataport);

  	/*Send our node info. to charmrun.
  	CommLock hasn't been initialized yet-- 
  	use non-locking version*/
  	ctrl_sendone_nolock("initnode",(const char *)&me,sizeof(me),NULL,0);
  
  	/*We get the other node addresses from a message sent
  	  back via the charmrun control port.*/
  	if (!skt_select1(Cmi_charmrun_fd,600*1000)) CmiAbort("Timeout waiting for nodetab!\n");
  	ChMessage_recv(Cmi_charmrun_fd,&nodetabmsg);
  }
  node_addresses_store(&nodetabmsg);
  ChMessage_free(&nodetabmsg);
}

#if CMK_NODE_QUEUE_AVAILABLE

/***********************************************************************
 * DeliverOutgoingNodeMessage()
 *
 * This function takes care of delivery of outgoing node messages from the
 * sender end. Broadcast messages are divided into sets of messages that 
 * are bound to the local node, and to remote nodes. For local
 * transmission, the messages are directly pushed into the recv
 * queues. For non-local transmission, the function DeliverViaNetwork()
 * is called
 ***********************************************************************/
void DeliverOutgoingNodeMessage(OutgoingMsg ogm)
{
  int i, rank, dst; OtherNode node;

  dst = ogm->dst;
  switch (dst) {
  case NODE_BROADCAST_ALL:
    CmiPushNode(CopyMsg(ogm->data,ogm->size));
    /*case-fallthrough (no break)-- deliver to all other processors*/
  case NODE_BROADCAST_OTHERS:
    for (i=0; i<_Cmi_numnodes; i++)
      if (i!=_Cmi_mynode)
	DeliverViaNetwork(ogm, nodes + i, DGRAM_NODEMESSAGE);
    GarbageCollectMsg(ogm);
    break;
  default:
    node = nodes+dst;
    rank=DGRAM_NODEMESSAGE;
    if (dst != _Cmi_mynode) {
      DeliverViaNetwork(ogm, node, rank);
      GarbageCollectMsg(ogm);
    } else {
      if (ogm->freemode == 'A') {
	CmiPushNode(CopyMsg(ogm->data,ogm->size));
	ogm->freemode = 'X';
      } else {
	CmiPushNode(ogm->data);
	FreeOutgoingMsg(ogm);
      }
    }
  }
}

#else

#define DeliverOutgoingNodeMessage(msg) DeliverOutgoingMessage(msg)

#endif

/***********************************************************************
 * DeliverOutgoingMessage()
 *
 * This function takes care of delivery of outgoing messages from the
 * sender end. Broadcast messages are divided into sets of messages that 
 * are bound to the local node, and to remote nodes. For local
 * transmission, the messages are directly pushed into the recv
 * queues. For non-local transmission, the function DeliverViaNetwork()
 * is called
 ***********************************************************************/
void DeliverOutgoingMessage(OutgoingMsg ogm)
{
  int i, rank, dst; OtherNode node;
  
  dst = ogm->dst;
  switch (dst) {
  case PE_BROADCAST_ALL:
    for (rank = 0; rank<_Cmi_mynodesize; rank++) {
      CmiPushPE(rank,CopyMsg(ogm->data,ogm->size));
    }
    for (i=0; i<_Cmi_numnodes; i++)
      if (i!=_Cmi_mynode)
	DeliverViaNetwork(ogm, nodes + i, DGRAM_BROADCAST);
    GarbageCollectMsg(ogm);
    break;
  case PE_BROADCAST_OTHERS:
    for (rank = 0; rank<_Cmi_mynodesize; rank++)
      if (rank + Cmi_nodestart != ogm->src) {
	CmiPushPE(rank,CopyMsg(ogm->data,ogm->size));
      }
    for (i = 0; i<_Cmi_numnodes; i++)
      if (i!=_Cmi_mynode)
	DeliverViaNetwork(ogm, nodes + i, DGRAM_BROADCAST);
    GarbageCollectMsg(ogm);
    break;
  default:
#ifndef CMK_OPTIMIZE
    if (dst<0 || dst>=CmiNumPes())
      CmiAbort("Send to out-of-bounds processor!");
#endif
    node = nodes_by_pe[dst];
    rank = dst - node->nodestart;
    if (node->nodestart != Cmi_nodestart) {
      DeliverViaNetwork(ogm, node, rank);
      GarbageCollectMsg(ogm);
    } else {
      if (ogm->freemode == 'A') {
	CmiPushPE(rank,CopyMsg(ogm->data,ogm->size));
	ogm->freemode = 'X';
      } else {
	CmiPushPE(rank, ogm->data);
	FreeOutgoingMsg(ogm);
      }
    }
  }
}


/******************************************************************************
 *
 * CmiGetNonLocal
 *
 * The design of this system is that the communication thread does all the
 * work, to eliminate as many locking issues as possible.  This is the only
 * part of the code that happens in the receiver-thread.
 *
 * This operation is fairly cheap, it might be worthwhile to inline
 * the code into CmiDeliverMsgs to reduce function call overhead.
 *
 *****************************************************************************/

#if CMK_NODE_QUEUE_AVAILABLE
char *CmiGetNonLocalNodeQ(void)
{
  char *result = 0;
  if(!PCQueueEmpty(CsvAccess(NodeState).NodeRecv)) {
    CmiLock(CsvAccess(NodeState).CmiNodeRecvLock);
    result = (char *) PCQueuePop(CsvAccess(NodeState).NodeRecv);
    CmiUnlock(CsvAccess(NodeState).CmiNodeRecvLock);
  }
  return result;
}
#endif

void *CmiGetNonLocal(void)
{
  CmiState cs = CmiGetState();
  CmiIdleLock_checkMessage(&cs->idle);
  return (void *) PCQueuePop(cs->recv);
}


/**
 * Set up and OutgoingMsg structure for this message.
 */
static OutgoingMsg PrepareOutgoing(CmiState cs,int pe,int size,int freemode,char *data) {
  OutgoingMsg ogm;
  MallocOutgoingMsg(ogm);
  CmiMsgHeaderSetLength(data, size);
  MACHSTATE2(2,"Preparing outgoing message for pe %d, size %d",pe,size);
  ogm->size = size;
  ogm->data = data;
  ogm->src = cs->pe;
  ogm->dst = pe;
  ogm->freemode = freemode;
  ogm->refcount = 0;
  return (CmiCommHandle)ogm;	
}

#if CMK_NODE_QUEUE_AVAILABLE

/******************************************************************************
 *
 * CmiGeneralNodeSend
 *
 * Description: This is a generic message sending routine. All the
 * converse message send functions are implemented in terms of this
 * function. (By setting appropriate flags (eg freemode) that tell
 * CmiGeneralSend() how exactly to handle the particular case of
 * message send)
 *
 *****************************************************************************/

CmiCommHandle CmiGeneralNodeSend(int node, int size, int freemode, char *data)
{
  
  CmiState cs = CmiGetState(); OutgoingMsg ogm;

  if (freemode == 'S') {
    char *copy = (char *)CmiAlloc(size);
    if (!copy)
      fprintf(stderr, "%d: Out of mem\n", _Cmi_mynode);
    memcpy(copy, data, size);
    data = copy; freemode = 'F';
  }

#if CMK_IMMEDIATE_MSG
    /* execute the immediate message right away */
  if (node == CmiMyNode() && CmiIsImmediate(data)) {
    CmiPushImmediateMsg(data);
    if (!_immRunning) CmiHandleImmediate();
    return;
  }
#endif
  ogm=PrepareOutgoing(cs,node,size,freemode,data);
  CmiCommLock();
  DeliverOutgoingNodeMessage(ogm);
  CmiCommUnlock();
  /* Check if any packets have arrived recently (preserves kernel network buffers). */
  CommunicationServer(0, 2);
  return (CmiCommHandle)ogm;
}

#endif


/******************************************************************************
 *
 * CmiGeneralSend
 *
 * Description: This is a generic message sending routine. All the
 * converse message send functions are implemented in terms of this
 * function. (By setting appropriate flags (eg freemode) that tell
 * CmiGeneralSend() how exactly to handle the particular case of
 * message send)
 *
 *****************************************************************************/

CmiCommHandle CmiGeneralSend(int pe, int size, int freemode, char *data)
{
  CmiState cs = CmiGetState(); OutgoingMsg ogm;

  if (freemode == 'S') {
    char *copy = (char *)CmiAlloc(size);
    if (!copy)
      fprintf(stderr, "%d: Out of mem\n", _Cmi_mynode);
    memcpy(copy, data, size);
    data = copy; freemode = 'F';
  }
  
  if (pe == cs->pe) {
#if ! CMK_SMP
    if (!_immRunning) /* CdsFifo_Enqueue, below, isn't SIGIO or thread safe.  
                      The SMP comm thread never gets here, because of the pe test. */
#endif
    {
#if CMK_IMMEDIATE_MSG
      /* execute the immediate message right away */
      /* but to avoid infinite recursive call, don't do this if _immRunning */
    if (CmiIsImmediate(data)) {
      CmiPushImmediateMsg(data);
      CmiHandleImmediate();
      return 0;
    }
#endif
    CdsFifo_Enqueue(cs->localqueue, data);
    if (freemode == 'A') {
      MallocOutgoingMsg(ogm);
      ogm->freemode = 'X';
      return ogm;
    } else return 0;
    }
  }

#if CMK_PERSISTENT_COMM
  if (phs) {
      CmiAssert(phsSize == 1);
      CmiSendPersistentMsg(*phs, pe, size, data);
      return NULL;
  }
#endif

  ogm=PrepareOutgoing(cs,pe,size,freemode,data);
  CmiCommLock();
  DeliverOutgoingMessage(ogm);
  CmiCommUnlock();
  /* Check if any packets have arrived recently (preserves kernel network buffers). */
  CommunicationServer(0, 2);
  return (CmiCommHandle)ogm;
}


void CmiSyncSendFn(int p, int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), 1);
  CmiGeneralSend(p,s,'S',m); 
}

CmiCommHandle CmiAsyncSendFn(int p, int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), 1);
  return CmiGeneralSend(p,s,'A',m); 
}

void CmiFreeSendFn(int p, int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), 1);
  CmiGeneralSend(p,s,'F',m); 
}

void CmiSyncBroadcastFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumPes()-1); 
  CmiGeneralSend(PE_BROADCAST_OTHERS,s,'S',m); 
}

CmiCommHandle CmiAsyncBroadcastFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumPes()-1); 
  return CmiGeneralSend(PE_BROADCAST_OTHERS,s,'A',m); 
}

void CmiFreeBroadcastFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumPes()-1);
  CmiGeneralSend(PE_BROADCAST_OTHERS,s,'F',m); 
}

void CmiSyncBroadcastAllFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumPes()); 
  CmiGeneralSend(PE_BROADCAST_ALL,s,'S',m); 
}

CmiCommHandle CmiAsyncBroadcastAllFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumPes()); 
  return CmiGeneralSend(PE_BROADCAST_ALL,s,'A',m); 
}

void CmiFreeBroadcastAllFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumPes()); 
  CmiGeneralSend(PE_BROADCAST_ALL,s,'F',m); 
}

#if CMK_NODE_QUEUE_AVAILABLE

void CmiSyncNodeSendFn(int p, int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), 1);
  CmiGeneralNodeSend(p,s,'S',m); 
}

CmiCommHandle CmiAsyncNodeSendFn(int p, int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), 1);
  return CmiGeneralNodeSend(p,s,'A',m); 
}

void CmiFreeNodeSendFn(int p, int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), 1);
  CmiGeneralNodeSend(p,s,'F',m); 
}

void CmiSyncNodeBroadcastFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumNodes()-1);
  CmiGeneralNodeSend(NODE_BROADCAST_OTHERS,s,'S',m); 
}

CmiCommHandle CmiAsyncNodeBroadcastFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumNodes()-1);
  return CmiGeneralNodeSend(NODE_BROADCAST_OTHERS,s,'A',m);
}

void CmiFreeNodeBroadcastFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumNodes()-1);
  CmiGeneralNodeSend(NODE_BROADCAST_OTHERS,s,'F',m); 
}

void CmiSyncNodeBroadcastAllFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumNodes());
  CmiGeneralNodeSend(NODE_BROADCAST_ALL,s,'S',m); 
}

CmiCommHandle CmiAsyncNodeBroadcastAllFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumNodes());
  return CmiGeneralNodeSend(NODE_BROADCAST_ALL,s,'A',m); 
}

void CmiFreeNodeBroadcastAllFn(int s, char *m)
{ 
  CQdCreate(CpvAccess(cQdState), CmiNumNodes());
  CmiGeneralNodeSend(NODE_BROADCAST_ALL,s,'F',m); 
}
#endif

/******************************************************************************
 *
 * Comm Handle manipulation.
 *
 *****************************************************************************/

int CmiAsyncMsgSent(CmiCommHandle handle)
{
  return (((OutgoingMsg)handle)->freemode == 'X');
}

void CmiReleaseCommHandle(CmiCommHandle handle)
{
  FreeOutgoingMsg(((OutgoingMsg)handle));
}

#if ! CMK_MULTICAST_LIST_USE_COMMON_CODE

/*****************************************************************************
 *
 * NET version List-Cast and Multicast Code
 *
 ****************************************************************************/
                                                                                
void CmiSyncListSendFn(int npes, int *pes, int len, char *msg)
{
  CmiError("ListSend not implemented.");
}
                                                                                
CmiCommHandle CmiAsyncListSendFn(int npes, int *pes, int len, char *msg)
{
  CmiError("ListSend not implemented.");
  return (CmiCommHandle) 0;
}
                                                                                
/* 
  because in all net versions, the message buffer after CmiSyncSendAndFree
  returns is not changed, we can use memory reference trick to avoid 
  memory copying here
*/
void CmiFreeListSendFn(int npes, int *pes, int len, char *msg)
{
  int i;
  for(i=0;i<npes;i++) {
    CmiReference(msg);
    CmiSyncSendAndFree(pes[i], len, msg);
  }
  CmiFree(msg);
}

#endif

/*
#if CMK_IMMEDIATE_MSG
void CmiProbeImmediateMsg()
{
  CommunicationServer(0, 0);
}
#endif
*/

/* Network progress function is used to poll the network when for
   messages. This flushes receive buffers on some implementations*/ 
void CmiMachineProgressImpl()
{
    CommunicationServer(0, 0);
}

/******************************************************************************
 *
 * Main code, Init, and Exit
 *
 *****************************************************************************/
extern void CthInit(char **argv);
extern void ConverseCommonInit(char **);

static char     **Cmi_argv;
static char     **Cmi_argvcopy;
static CmiStartFn Cmi_startfn;   /* The start function */
static int        Cmi_usrsched;  /* Continue after start function finishes? */

static void ConverseRunPE(int everReturn)
{
  CmiIdleState *s=CmiNotifyGetState();
  CmiState cs;
  char** CmiMyArgv;
  CmiNodeAllBarrier();
  cs = CmiGetState();
  CpvInitialize(void *,CmiLocalQueue);
  CpvAccess(CmiLocalQueue) = cs->localqueue;

  /* all non 0 rank use the copied one while rank 0 will modify the actual argv */
  if (CmiMyRank())
    CmiMyArgv = CmiCopyArgs(Cmi_argvcopy);
  else
    CmiMyArgv = Cmi_argv;
  CthInit(CmiMyArgv);

#if CMK_USE_GM
  CmiCheckGmStatus();
#endif

  ConverseCommonInit(CmiMyArgv);

  /* initialize the network progress counter*/
  /* Network progress function is used to poll the network when for
     messages. This flushes receive buffers on some  implementations*/ 
  CpvInitialize(int , networkProgressCount);
  CpvAccess(networkProgressCount) = 0;

  /* better to show the status here */
  if (Cmi_netpoll == 1 && CmiMyPe() == 0)
    CmiPrintf("Charm++: scheduler running in netpoll mode.\n");
  
#if CMK_USE_GM
  if (Cmi_charmrun_fd != -1)
#endif
  {
  CcdCallOnConditionKeep(CcdPROCESSOR_BEGIN_IDLE,
      (CcdVoidFn) CmiNotifyBeginIdle, (void *) s);
  CcdCallOnConditionKeep(CcdPROCESSOR_STILL_IDLE,
      (CcdVoidFn) CmiNotifyStillIdle, (void *) s);
  }
#if CMK_SHARED_VARS_UNAVAILABLE
  if (Cmi_netpoll) /*Repeatedly call CommServer*/
    CcdCallOnConditionKeep(CcdPERIODIC, 
        (CcdVoidFn) CommunicationPeriodic, NULL);
  else /*Only need this for retransmits*/
    CcdCallOnConditionKeep(CcdPERIODIC_10ms, 
        (CcdVoidFn) CommunicationPeriodic, NULL);
#endif

  if (CmiMyRank()==0 && Cmi_charmrun_fd!=-1) {
    CcdCallOnConditionKeep(CcdPERIODIC_10ms, (CcdVoidFn) CmiStdoutFlush, NULL);
#if CMK_SHARED_VARS_UNAVAILABLE
    if (!Cmi_asyncio) {
    /* gm cannot live with setitimer */
    CcdCallFnAfter((CcdVoidFn)pingCharmrunPeriodic,NULL,1000);
    }
    else {
    /*Occasionally ping charmrun, to test if it's dead*/
    struct itimerval i;
    CmiSignal(SIGALRM, 0, 0, pingCharmrun);
    i.it_interval.tv_sec = 1;
    i.it_interval.tv_usec = 0;
    i.it_value.tv_sec = 1;
    i.it_value.tv_usec = 0;
    setitimer(ITIMER_REAL, &i, NULL);
    }

#if ! CMK_USE_GM && ! CMK_USE_TCP
    /*Occasionally check for retransmissions, outgoing acks, etc.*/
    /*no need for GM case */
    CcdCallFnAfter((CcdVoidFn)CommunicationsClockCaller,NULL,Cmi_comm_clock_delay);
#endif
#endif
    
    /*Initialize the clock*/
    Cmi_clock=GetClock();
  }
  
#ifdef CMK_RANDOMLY_CORRUPT_MESSAGES
  srand((int)(1024.0*CmiWallTimer()));
  if (CmiMyPe()==0)
    CmiPrintf("Charm++: Machine layer will randomly corrupt every %d'th message (rand %d)\n",
    	CMK_RANDOMLY_CORRUPT_MESSAGES,rand());
#endif

  /* communication thread */
  if (CmiMyRank() == CmiMyNodeSize()) {
    Cmi_startfn(CmiGetArgc(CmiMyArgv), CmiMyArgv);
    if (Cmi_charmrun_fd!=-1)
          while (1) CommunicationServer(5, 0);
  }
  else
  if (!everReturn) {
    Cmi_startfn(CmiGetArgc(CmiMyArgv), CmiMyArgv);
    /* turn on immediate messages only now
     node barrier previously should take care of the node synchronization */
    _immediateReady = 1;
    if (Cmi_usrsched==0) CsdScheduler(-1);
    ConverseExit();
  }
}

void ConverseExit(void)
{
  machine_initiated_shutdown=1;
  if (CmiMyRank()==0) {
    if(Cmi_print_stats)
      printNetStatistics();
    log_done();
    CmiMachineExit();
    CmiStdoutFlush();
    ConverseCommonExit();
    if (Cmi_charmrun_fd==-1) exit(0); /*Standalone version-- just leave*/
  }
  if (Cmi_charmrun_fd!=-1) {
  	ctrl_sendone_locking("ending",NULL,0,NULL,0); /* this causes charmrun to go away */
#if CMK_SHARED_VARS_UNAVAILABLE
 	Cmi_check_delay = 1.0;		/* speed up checking of charmrun */
 	while (1) CommunicationServer(500, 2);
#endif
  }
/*Comm. thread will kill us.*/
  while (1) CmiYield();
}

static void set_signals(void)
{
  if(!Cmi_truecrash) {
    signal(SIGSEGV, KillOnAllSigs);
    signal(SIGFPE, KillOnAllSigs);
    signal(SIGILL, KillOnAllSigs);
    signal(SIGINT, KillOnAllSigs);
    signal(SIGTERM, KillOnAllSigs);
    signal(SIGABRT, KillOnAllSigs);
#   if !defined(_WIN32) || defined(__CYGWIN__) /*UNIX-only signals*/
    signal(SIGQUIT, KillOnAllSigs);
    signal(SIGBUS, KillOnAllSigs);
#     if CMK_HANDLE_SIGUSR
    signal(SIGUSR1, HandleUserSignals);
    signal(SIGUSR2, HandleUserSignals);
#     endif
#   endif /*UNIX*/
  }
}

/*Socket idle function to use before addresses have been
  obtained.  During the real program, we idle with CmiYield.
*/
static void obtain_idleFn(void) {sleep(0);}

static int net_default_skt_abort(int code,const char *msg)
{
  fprintf(stderr,"Fatal socket error: code %d-- %s\n",code,msg);
  machine_exit(1);
  return -1;
}

void ConverseInit(int argc, char **argv, CmiStartFn fn, int usc, int everReturn)
{
#if MACHINE_DEBUG
  debugLog=NULL;
#endif
#if CMK_USE_HP_MAIN_FIX
#if FOR_CPLUS
  _main(argc,argv);
#endif
#endif
  Cmi_argvcopy = CmiCopyArgs(argv);
  Cmi_argv = argv; Cmi_startfn = fn; Cmi_usrsched = usc;
  Cmi_netpoll = 0;
#if CMK_NETPOLL
  Cmi_netpoll = 1;
#endif
#if CMK_WHEN_PROCESSOR_IDLE_USLEEP
  Cmi_idlepoll = 0;
#else
  Cmi_idlepoll = 1;
#endif
  Cmi_truecrash = 0;
  if (CmiGetArgFlagDesc(argv,"+truecrash","Do not install signal handlers") ||
      CmiGetArgFlagDesc(argv,"++debug",NULL /*meaning: don't show this*/)) Cmi_truecrash = 1;
    /* netpoll disable signal */
  if (CmiGetArgFlagDesc(argv,"+netpoll","Do not use SIGIO--poll instead")) Cmi_netpoll = 1;
    /* idlepoll use poll instead if sleep when idle */
  if (CmiGetArgFlagDesc(argv,"+idlepoll","Do not sleep when idle")) Cmi_idlepoll = 1;
    /* idlesleep use sleep instead if busywait when idle */
  if (CmiGetArgFlagDesc(argv,"+idlesleep","Make sleep calls when idle")) Cmi_idlepoll = 0;
  Cmi_syncprint = CmiGetArgFlagDesc(argv,"+syncprint", "Flush each CmiPrintf to the terminal");

  Cmi_asyncio = 1;
#if CMK_ASYNC_NOT_NEEDED
  Cmi_asyncio = 0;
#endif
  if (CmiGetArgFlagDesc(argv,"+asyncio","Use async IO")) Cmi_asyncio = 1;
  if (CmiGetArgFlagDesc(argv,"+asynciooff","Don not use async IO")) Cmi_asyncio = 0;

  MACHSTATE2(5,"Init: (netpoll=%d), (idlepoll=%d)",Cmi_netpoll,Cmi_idlepoll);

  skt_init();
  /* use special abort handler instead of default_skt_abort to 
     prevent exit trapped by atexit_check() due to the exit() call  */
  skt_set_abort(net_default_skt_abort);
  atexit(machine_atexit_check);
  parse_netstart();
  extract_args(argv);
  log_init();
  Cmi_scanf_mutex = CmiCreateLock();

#if MACHINE_DEBUG_LOG
  {
    char ln[200];
    sprintf(ln,"debugLog.%d",_Cmi_mynode);
    debugLog=fopen(ln,"w");
  }
#endif

  skt_set_idle(obtain_idleFn);
  if (!skt_ip_match(Cmi_charmrun_IP,_skt_invalid_ip)) {
  	set_signals();
#if CMK_USE_TCP
  	dataskt=skt_server(&dataport);
#elif !CMK_USE_GM
  	dataskt=skt_datagram(&dataport, Cmi_os_buffer_size);
#else
        dataskt=-1;
#endif
	MACHSTATE2(5,"skt_connect at dataskt:%d Cmi_charmrun_port:%d",dataskt, Cmi_charmrun_port);
  	Cmi_charmrun_fd = skt_connect(Cmi_charmrun_IP, Cmi_charmrun_port, 1800);
	MACHSTATE1(5,"Opened data socket at port %d",dataport);
	CmiStdoutInit();
  } else {/*Standalone operation*/
  	printf("Charm++: standalone mode (not using charmrun)\n");
  	dataskt=-1;
  	Cmi_charmrun_fd=-1;
  }

  CmiMachineInit(argv);

  node_addresses_obtain(argv);

#if CMK_USE_TCP
  open_tcp_sockets();
#endif

  skt_set_idle(CmiYield);
  Cmi_check_delay = 2.0+0.5*_Cmi_numnodes;
  if (Cmi_charmrun_fd==-1) /*Don't bother with check in standalone mode*/
	Cmi_check_delay=1.0e30;

  CsvInitialize(CmiNodeState, NodeState);
  CmiNodeStateInit(&CsvAccess(NodeState));

  /* Network progress function is used to poll the network when for
     messages. This flushes receive buffers on some  implementations*/ 
  networkProgressPeriod = 0;  
  CmiGetArgInt(argv, "+networkProgressPeriod", &networkProgressPeriod);

  CmiStartThreads(argv);
  ConverseRunPE(everReturn);
}


#if CMK_PERSISTENT_COMM

int persistentSendMsgHandlerIdx;

static void sendPerMsgHandler(char *msg)
{
  int msgSize;
  void *destAddr, *destSizeAddr;

  msgSize = CmiMsgHeaderGetLength(msg);
  msgSize -= 2*sizeof(void *);
  destAddr = *(void **)(msg + msgSize);
  destSizeAddr = *(void **)(msg + msgSize + sizeof(void*));
/*CmiPrintf("msgSize:%d destAddr:%p, destSizeAddr:%p\n", msgSize, destAddr, destSizeAddr);*/
  CmiSetHandler(msg, CmiGetXHandler(msg));
  *((int *)destSizeAddr) = msgSize;
  memcpy(destAddr, msg, msgSize);
}

void CmiSendPersistentMsg(PersistentHandle h, int destPE, int size, void *m)
{
  CmiAssert(h!=NULL);
  PersistentSendsTable *slot = (PersistentSendsTable *)h;
  CmiAssert(slot->used == 1);
  CmiAssert(slot->destPE == destPE);
  if (size > slot->sizeMax) {
    CmiPrintf("size: %d sizeMax: %d\n", size, slot->sizeMax);
    CmiAbort("Abort: Invalid size\n");
  }

/*CmiPrintf("[%d] CmiSendPersistentMsg h=%p hdl=%d destpe=%d destAddress=%p size=%d\n", CmiMyPe(), *phs, CmiGetHandler(m), slot->destPE, slot->destAddress, size);*/

  if (slot->destAddress[0]) {
    int newsize = size + sizeof(void *)*2;
    char *newmsg = (char*)CmiAlloc(newsize);
    memcpy(newmsg, m, size);
    memcpy(newmsg+size, &slot->destAddress[0], sizeof(void *));
    memcpy(newmsg+size+sizeof(void*), &slot->destSizeAddress[0], sizeof(void *));
    CmiFree(m);
    CmiMsgHeaderSetLength(newmsg, size + sizeof(void *)*2);
    CmiSetXHandler(newmsg, CmiGetHandler(newmsg));
    CmiSetHandler(newmsg, persistentSendMsgHandlerIdx);
    phs = NULL; phsSize = 0;
    CmiSyncSendAndFree(slot->destPE, newsize, newmsg);
  }
  else {
#if 1
    if (slot->messageBuf != NULL) {
      CmiPrintf("Unexpected message in buffer on %d\n", CmiMyPe());
      CmiAbort("");
    }
    slot->messageBuf = m;
    slot->messageSize = size;
#else
    /* normal send */
    PersistentHandle  *phs_tmp = phs;
    int phsSize_tmp = phsSize;
    phs = NULL; phsSize = 0;
    CmiPrintf("[%d]Slot sending message directly\n", CmiMyPe());
    CmiSyncSendAndFree(slot->destPE, size, m);
    phs = phs_tmp; phsSize = phsSize_tmp;
#endif
  }
}

void CmiSyncSendPersistent(int destPE, int size, char *msg, PersistentHandle h)
{
  char *dupmsg = (char *) CmiAlloc(size);
  memcpy(dupmsg, msg, size);

  /*  CmiPrintf("Setting root to %d\n", 0); */
  if (CmiMyPe()==destPE) {
    CQdCreate(CpvAccess(cQdState), 1);
    CdsFifo_Enqueue(CpvAccess(CmiLocalQueue),dupmsg);
  }
  else
    CmiSendPersistentMsg(h, destPE, size, dupmsg);
}

/* called in PumpMsgs */
int PumpPersistent()
{
  PersistentReceivesTable *slot = persistentReceivesTableHead;
  int status = 0;
  while (slot) {
    unsigned int size = *(slot->recvSizePtr[0]);
    if (size > 0)
    {
      char *msg = slot->messagePtr[0];
/*CmiPrintf("size: %d msg:%p %p\n", size, msg, slot->messagePtr);*/

#if 0
      void *dupmsg;
      dupmsg = CmiAlloc(size);
      
      _MEMCHECK(dupmsg);
      memcpy(dupmsg, msg, size);
      msg = dupmsg;
#else
      /* return messagePtr directly and user MUST make sure not to delete it. */
      /*CmiPrintf("[%d] %p size:%d rank:%d root:%d\n", CmiMyPe(), msg, size, CMI_DEST_RANK(msg), CMI_BROADCAST_ROOT(msg));*/

      CmiReference(msg);
#endif

      CmiPushPE(CMI_DEST_RANK(msg), msg);
#if CMK_BROADCAST_SPANNING_TREE
      if (CMI_BROADCAST_ROOT(msg))
          SendSpanningChildren(size, msg);
#endif
      *(slot->recvSizePtr[0]) = 0;
      status = 1;
    }
    slot = slot->next;
  }
  return status;
}

void *PerAlloc(int size)
{
  return CmiAlloc(size);
}
                                                                                
void PerFree(char *msg)
{
    CmiFree(msg);
}

void persist_machine_init() 
{
  persistentSendMsgHandlerIdx =
       CmiRegisterHandler((CmiHandler)sendPerMsgHandler);
}

#endif

/*@}*/
