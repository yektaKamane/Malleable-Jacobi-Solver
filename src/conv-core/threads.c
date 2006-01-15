/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

 /**************************************************************************
 *
 * typedef CthThread
 *
 *   - a first-class thread object.
 *
 * CthThread CthSelf()
 *
 *   - returns the current thread.
 *
 * void CthResume(CthThread t)
 *
 *   - Immediately transfers control to thread t.  Note: normally, the user
 *     of a thread package wouldn't explicitly choose which thread to transfer
 *     to.  Instead, the user would rely upon a "scheduler" to choose the
 *     next thread.  Therefore, this routine is primarily intended for people
 *     who are implementing schedulers, not for end-users.  End-users should
 *     probably call CthSuspend or CthAwaken (see below).
 *
 * CthThread CthCreate(CthVoidFn fn, void *arg, int size)
 *
 *   - Creates a new thread object.  The thread is not given control yet.
 *     The thread is not passed to the scheduler.  When (and if) the thread
 *     eventually receives control, it will begin executing the specified 
 *     function 'fn' with the specified argument.  The 'size' parameter
 *     specifies the stack size, 0 means use the default size.
 *
 * void CthFree(CthThread t)
 *
 *   - Frees thread t.  You may free the currently-executing thread, although
 *     the free will actually be postponed until the thread suspends.
 *
 *
 * In addition to the routines above, the threads package assumes that there
 * will be a "scheduler" of some sort, whose job is to select which threads
 * to execute.  The threads package does not provide a scheduler (although
 * converse may provide one or more schedulers separately).  However, for
 * standardization reasons, it does define an interface to which all schedulers
 * can comply.  A scheduler consists of a pair of functions:
 *
 *   - An awaken-function.  The awaken-function is called to
 *     notify the scheduler that a particular thread needs the CPU.  The
 *     scheduler is expected to respond to this by inserting the thread
 *     into a ready-pool of some sort.
 *
 *   - A choose-next function.  The choose-next function is called to
 *     to ask the scheduler which thread to execute next.
 *
 * The interface to the scheduler is formalized in the following functions:
 *
 * void CthSchedInit()
 *
 *   - you must call this before any of the following functions will work.
 *
 * void CthSuspend()
 *
 *   - The thread calls this function, which in turn calls the scheduler's
 *     choose-next function.  It then resumes whatever thread is returned
 *     by the choose-next function.
 *
 * void CthAwaken(CthThread t)
 *
 *   - The thread-package user calls this function, which in turn calls the
 *     scheduler's awaken-function to awaken thread t.  This probably causes
 *     the thread t to be inserted in the ready-pool.
 *
 * void CthSetStrategy(CthThread t, CthAwkFn awakenfn, CthThFn choosefn)
 *
 *     This specifies the scheduling functions to be used for thread 't'.
 *     The scheduling functions must have the following prototypes:
 *
 *          void awakenfn(CthThread t);
 *          CthThread choosefn();
 *
 *     These functions must be provided on a per-thread basis.  (Eg, if you
 *     CthAwaken a thread X, then X's awakefn will be called.  If a thread Y
 *     calls CthSuspend, then Y's choosefn will be called to pick the next
 *     thread.)  Of course, you may use the same functions for all threads
 *     (the common case), but the specification on a per-thread basis gives
 *     you maximum flexibility in controlling scheduling.
 *
 *     See also: common code, CthSetStrategyDefault.
 *
 * void CthYield()
 *
 *   - simply executes { CthAwaken(CthSelf()); CthSuspend(); }.  This
 *     combination gives up control temporarily, but ensures that control
 *     will eventually return.
 *
 *
 * Note: there are several possible ways to implement threads.   
 *****************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if CMK_MEMORY_PROTECTABLE
#include <malloc.h> /*<- for memalign*/
#endif
#if CMK_VERSION_BLUEGENE
#include "rts.h"	/*<- for rts_memory_alias */
#endif

#include "converse.h"
#include "qt.h"

#include "conv-trace.h"
#include <sys/types.h>

#ifndef CMK_STACKSIZE_DEFAULT
#define CMK_STACKSIZE_DEFAULT 32768
#endif

#if ! CMK_THREADS_BUILD_DEFAULT
#undef CMK_THREADS_USE_JCONTEXT
#undef CMK_THREADS_USE_CONTEXT
#undef CMK_THREADS_ARE_WIN32_FIBERS
#undef CMK_THREADS_USE_PTHREADS

#if CMK_THREADS_BUILD_CONTEXT
#define CMK_THREADS_USE_CONTEXT       1
#elif CMK_THREADS_BUILD_JCONTEXT
#define CMK_THREADS_USE_JCONTEXT       1
#elif  CMK_THREADS_BUILD_FIBERS
#define CMK_THREADS_ARE_WIN32_FIBERS  1
#elif  CMK_THREADS_BUILD_PTHREADS
#define CMK_THREADS_USE_PTHREADS      1
#elif  CMK_THREADS_BUILD_STACKCOPY
#define CMK_THREADS_USE_STACKCOPY      1
#endif

#endif

/**************************** Shared Base Thread Class ***********************/
/*
	FAULT_EVAC
	Moved the cmicore converse header from CthThreadBase to CthThreadToken.
	The CthThreadToken gets enqueued in the converse queue instead of the
	CthThread. This allows the thread to be moved out of a processor even
	if there is an awaken call for it enqueued in the scheduler

*/

typedef struct CthThreadBase
{
  CthThreadToken *token; /* token that shall be enqueued into the ready queue*/
  int scheduled;	 /* has this thread been added to the ready queue ? */
 
  CmiObjId   tid;        /* globally unique tid */
  CthAwkFn   awakenfn;   /* Insert this thread into the ready queue */
  CthThFn    choosefn;   /* Return the next ready thread */
  CthThread  next; /* Next active thread */
  int        suspendable; /* Can this thread be blocked */
  int        exiting;    /* Is this thread finished */

  char      *data;       /* thread private data */
  int        datasize;   /* size of thread-private data, in bytes */

  int isMigratable; /* thread is migratable (isomalloc or aliased stack) */
  int aliasStackHandle; /* handle for aliased stack */
  void      *stack; /*Pointer to thread stack*/
  int        stacksize; /*Size of thread stack (bytes)*/
  struct CthThreadListener *listener; /* pointer to the first of the listeners */
} CthThreadBase;

/* By default, there are no flags */
static int CmiThreadIs_flag=0;

int CmiThreadIs(int flag)
{
	return (CmiThreadIs_flag&flag)==flag;
}

/*Macros to convert between base and specific thread types*/
#define B(t) ((CthThreadBase *)(t))
#define S(t) ((CthThread)(t))


CthThreadToken *CthGetToken(CthThread t){
	return B(t)->token;
}

/*********************** Stack Aliasing *********************
  Stack aliasing: instead of consuming virtual address space
  with isomalloc, map all stacks to the same virtual addresses
  (like stack copying), but use the VM hardware to make thread
  swaps fast.  This implementation uses *files* to store the stacks,
  and mmap to drop the stack data into the right address.  Despite
  the presence of files, at least under Linux with local disks context
  switches are still fast; the mmap overhead is less than 5us per
  switch, even for multi-megabyte stacks.
  
  WARNING: Does NOT work in SMP mode, because all the processors
      on a node will try to map their stacks to the same place.
  WARNING: Does NOT work if switching directly from one migrateable
      thread to another, because it blows away the stack you're currently
      running on.  The usual Charm approach of switching to the 
      (non-migratable) scheduler thread on context switch works fine.
*/
#ifndef CMK_THREADS_ALIAS_STACK
#define CMK_THREADS_ALIAS_STACK 0
#endif

/** Address to map all migratable thread stacks to. */
#define CMK_THREADS_ALIAS_LOCATION ((void *)0xb0000000)

#if CMK_THREADS_ALIAS_STACK
#include <stdlib.h> /* for mkstemp */
#include <sys/mman.h> /* for mmap */
#include <errno.h> /* for perror */
#include <unistd.h> /* for unlink, lseek, close */

/** Create an aliasable area of this size.  Returns alias handle. */
int CthAliasCreate(int stackSize)
{
  /* Make a file to contain thread stack */
  char tmpName[1024];
  char lastByte=0;
  int fd;
  sprintf(tmpName,"/tmp/charmThreadStackXXXXXX");
  fd=mkstemp(tmpName);
  if (fd==-1) CmiAbort("threads.c> Cannot create /tmp file to contain thread stack");
  unlink(tmpName); /* delete file when it gets closed */
  
  /* Make file big enough for stack, by writing one byte at end */
  lseek(fd,stackSize-sizeof(lastByte),SEEK_SET);
  write(fd,&lastByte,sizeof(lastByte));
  
  return fd;
}

void CthAliasFree(int fd) {
  close(fd);
}

#endif

/**
 CthAliasEnable brings this thread's stack into memory.
 You must call it before accessing the thread stack, 
 for example, before running, packing, or unpacking the stack data.
*/
#if CMK_THREADS_ALIAS_STACK
CthThreadBase *_curMappedStack=0;
void CthAliasEnable(CthThreadBase *t) {
	void *s;
	int flags=MAP_FIXED|MAP_SHARED; /* Posix flags */
	if (!t->isMigratable) return;
	if (t==_curMappedStack) return; /* don't re-map */
	_curMappedStack=t;
	if (0) printf("Mmapping in thread %p from runtime stack %p\n",t,&s);
	
#if CMK_VERSION_BLUEGENE
        /* Blue Gene/L does not have mmap */
        /* So this depends on a hack in CNK for syscall rts_memory_alias to
           map stack pointer to a fix address */
        {
        register char *dest=(char *)0x70000000; /* point stack at high memory, not user memory */
        int alias_slot;
        rts_memory_alias(t->stack,dest,t->stacksize,&alias_slot);
        }
#else
	/* Linux mmap flag MAP_POPULATE, to pre-fault in all the pages,
	   only seems to slow down overall performance. */
	/* Linux mmap flag MAP_GROWSDOWN is rejected at runtime under 2.4.25 */
	s=mmap(t->stack,t->stacksize,
		PROT_READ|PROT_WRITE|PROT_EXEC, /* exec for gcc nested function thunks */
		flags, t->aliasStackHandle,0);
	if (s!=t->stack) {
		perror("threads.c CthAliasEnable mmap");
		CmiAbort("threads.c CthAliasEnable mmap failed");
	}
#endif
}
#else
#define CthAliasEnable(t) /* empty */
#endif





/*********** Thread-local storage *********/

CthCpvStatic(CthThread,  CthCurrent); /*Current thread*/
CthCpvDeclare(char *,    CthData); /*Current thread's private data (externally visible)*/
CthCpvStatic(int,        CthDatasize);

void CthSetThreadID(CthThread th, int a, int b, int c)
{
  B(th)->tid.id[0] = a;
  B(th)->tid.id[1] = b;
  B(th)->tid.id[2] = c;
}

/* possible hack? CW */
CmiObjId *CthGetThreadID(CthThread th)
{
  return &(B(th)->tid);
}

char *CthGetData(CthThread t) { return B(t)->data; }

/* Ensure this thread has at least enough 
room for all the thread-local variables 
initialized so far on this processor.
*/
static void CthFixData(CthThread t)
{
  int newsize = CthCpvAccess(CthDatasize);
  int oldsize = B(t)->datasize;
  if (oldsize < newsize) {
    newsize = 2*newsize;
    B(t)->datasize = newsize;
    /* Note: realloc(NULL,size) is equivalent to malloc(size) */
    B(t)->data = (char *)realloc(B(t)->data, newsize);
    memset(B(t)->data+oldsize, 0, newsize-oldsize);
  }
}

/**
Allocate another size bytes of thread-local storage,
and return the offset into the thread storage buffer.
 */
int CthRegister(int size)
{
  int datasize=CthCpvAccess(CthDatasize);
  CthThreadBase *th=(CthThreadBase *)CthCpvAccess(CthCurrent);
  int result, align = 1;
  while (size>align) align<<=1;
  datasize = (datasize+align-1) & ~(align-1);
  result = datasize;
  datasize += size;
  CthCpvAccess(CthDatasize) = datasize;
  CthFixData(S(th)); /*Make the current thread have this much storage*/
  CthCpvAccess(CthData) = th->data;
  return result;
}

/**
Make sure we have room to store up to at least maxOffset
bytes of thread-local storage.
*/
void CthRegistered(int maxOffset) {
  if (CthCpvAccess(CthDatasize)<maxOffset) {
    CthThreadBase *th=(CthThreadBase *)CthCpvAccess(CthCurrent);
    CthCpvAccess(CthDatasize) = maxOffset;
    CthFixData(S(th)); /*Make the current thread have this much storage*/
    CthCpvAccess(CthData) = th->data;
  }
}

/*********** Creation and Deletion **********/
CthCpvStatic(int, _defaultStackSize);

static void CthThreadBaseInit(CthThreadBase *th)
{
  static int serialno = 1;
#if CMK_LINUX_PTHREAD_HACK
  /*HACK for LinuxThreads: to support our user-level threads
    library, we use a slightly modified version of libpthread.a
    with user-level threads support enabled via these flags.
  */
  extern int __pthread_find_self_with_pid;
  extern int __pthread_nonstandard_stacks;
  __pthread_find_self_with_pid=1;
  __pthread_nonstandard_stacks=1;
#endif
	th->token = (CthThreadToken *)malloc(sizeof(CthThreadToken));
	th->token->thread = S(th);
	th->scheduled = 0;

  th->awakenfn = 0;
  th->choosefn = 0;
  th->next=0;
  th->suspendable = 1;
  th->exiting = 0;

  th->data=0;
  th->datasize=0;
  CthFixData(S(th));

  CthSetStrategyDefault(S(th));

  th->isMigratable=0;
  th->aliasStackHandle=0;
  th->stack=NULL;
  th->stacksize=0;

  th->tid.id[0] = CmiMyPe();
  th->tid.id[1] = serialno++;
  th->tid.id[2] = 0;

  th->listener = NULL;
}

static void *CthAllocateStack(CthThreadBase *th,int *stackSize,int useMigratable)
{
  void *ret=NULL;
  if (*stackSize==0) *stackSize=CthCpvAccess(_defaultStackSize);
  th->stacksize=*stackSize;
  if (!useMigratable) {
    ret=malloc(*stackSize); 
  } else {
    th->isMigratable=1;
#if CMK_THREADS_ALIAS_STACK
    th->aliasStackHandle=CthAliasCreate(*stackSize);
    ret=CMK_THREADS_ALIAS_LOCATION;
#else /* isomalloc */
    ret=CmiIsomalloc(*stackSize);
#endif
  }
  _MEMCHECK(ret);
  th->stack=ret;
  return ret;
}
static void CthThreadBaseFree(CthThreadBase *th)
{
	struct CthThreadListener *l,*lnext;
	/*
	 * remove the token if it is not queued in the converse scheduler		
	 */
	if(th->scheduled == 0){
		free(th->token);
	}else{
		th->token->thread = NULL;
	}
	/* Call the free function pointer on all the listeners on
			this thread and also delete the thread listener objects
	*/
	for(l=th->listener;l!=NULL;l=lnext){
			lnext=l->next;
			l->next=0;
			l->free(l);
	}
  free(th->data);
  if (th->isMigratable) {
#if CMK_THREADS_ALIAS_STACK
	  CthAliasFree(th->aliasStackHandle);
#else /* isomalloc */
	  CmiIsomallocFree(th->stack);
#endif
  } 
  else if (th->stack!=NULL) {
	  free(th->stack);
  }
  th->stack=NULL;
}

CpvDeclare(int, _numSwitches); /*Context switch count*/

static void CthBaseInit(char **argv)
{
  CpvInitialize(int, _numSwitches);
  CpvAccess(_numSwitches) = 0;

  CthCpvInitialize(int,  _defaultStackSize);
  CthCpvAccess(_defaultStackSize)=CMK_STACKSIZE_DEFAULT;
  CmiGetArgIntDesc(argv,"+stacksize",&CthCpvAccess(_defaultStackSize),
  	"Default user-level thread stack size");  
  
  CthCpvInitialize(CthThread,  CthCurrent);
  CthCpvInitialize(char *, CthData);
  CthCpvInitialize(int,        CthDatasize);
  
  CthCpvAccess(CthData)=0;
  CthCpvAccess(CthDatasize)=0;
}

int CthImplemented() { return 1; } 

CthThread CthSelf()
{
  return CthCpvAccess(CthCurrent);
}

void CthPupBase(pup_er p,CthThreadBase *t,int useMigratable)
{
#ifndef CMK_OPTIMIZE
	if ((CthThread)t==CthCpvAccess(CthCurrent))
		CmiAbort("CthPupBase: Cannot pack running thread!");
#endif
	/*recreate the token when unpacking and set the thread
	  pointer in it to NULL if the thread is being deleted*/
	if(pup_isUnpacking(p)){
		t->token = (CthThreadToken *)malloc(sizeof(CthThreadToken));
		t->token->thread = S(t);
		t->scheduled = 0;
	}
	/*Really need a pup_functionPtr here:*/
	pup_bytes(p,&t->awakenfn,sizeof(t->awakenfn));
	pup_bytes(p,&t->choosefn,sizeof(t->choosefn));
	pup_bytes(p,&t->next,sizeof(t->next));
	pup_int(p,&t->suspendable);
	pup_int(p,&t->datasize);
	if (pup_isUnpacking(p)) { 
		t->data = (char *) malloc(t->datasize);_MEMCHECK(t->data);
	}
	pup_bytes(p,(void *)t->data,t->datasize);
	pup_int(p,&t->isMigratable);
	pup_int(p,&t->stacksize);
	if (t->isMigratable) {
#if CMK_THREADS_ALIAS_STACK
		if (pup_isUnpacking(p)) { 
			CthAllocateStack(t,&t->stacksize,1);
		}
		CthAliasEnable(t);
		pup_bytes(p,t->stack,t->stacksize);
#else /* isomalloc */
		CmiIsomallocPup(p,&t->stack);
#endif
	} 
	else {
		if (useMigratable)
			CmiAbort("You must use CthCreateMigratable to use CthPup!\n");
		/*Pup the stack pointer as raw bytes*/
		pup_bytes(p,&t->stack,sizeof(t->stack));
	}
	if (pup_isUnpacking(p)) { 
		/* FIXME:  restore thread listener */
		t->listener = NULL;
	}
}

static void CthThreadFinished(CthThread t)
{
	B(t)->exiting=1;
	CthSuspend();
}

/************ Scheduler Interface **********/

void CthSetSuspendable(CthThread t, int val) { B(t)->suspendable = val; }
int CthIsSuspendable(CthThread t) { return B(t)->suspendable; }

void CthSetNext(CthThread t, CthThread v) { B(t)->next = v; }
CthThread CthGetNext(CthThread t) { return B(t)->next; }

static void CthNoStrategy(void)
{
  CmiAbort("Called CthAwaken or CthSuspend before calling CthSetStrategy.\n");
}

void CthSetStrategy(CthThread t, CthAwkFn awkfn, CthThFn chsfn)
{
  B(t)->awakenfn = awkfn;
  B(t)->choosefn = chsfn;
}

static void CthBaseResume(CthThread t)
{
	
	struct CthThreadListener *l;
	for(l=B(t)->listener;l!=NULL;l=l->next){
			l->resume(l);
	}
  CpvAccess(_numSwitches)++;
  CthFixData(t); /*Thread-local storage may have changed in other thread.*/
  CthCpvAccess(CthCurrent) = t;
  CthCpvAccess(CthData) = B(t)->data;
  CthAliasEnable(B(t));
}

/**
  switch the thread to t
*/
void CthSwitchThread(CthThread t)
{
  CthBaseResume(t);
}

/*
Suspend: finds the next thread to execute, and resumes it
*/
void CthSuspend(void)
{
  CthThread next;
	struct CthThreadListener *l;
  CthThreadBase *cur=B(CthCpvAccess(CthCurrent));

  if (cur->suspendable == 0)
    CmiAbort("Fatal Error> trying to suspend a non-suspendable thread!\n");
  
	/*
		Call the suspend function on listeners
	*/
	for(l=cur->listener;l!=NULL;l=l->next){
			l->suspend(l);
	}
  if (cur->choosefn == 0) CthNoStrategy();
  next = cur->choosefn();
	cur->scheduled=0;
#ifndef CMK_OPTIMIZE
#if !CMK_TRACE_IN_CHARM
  if(CpvAccess(traceOn))
    traceSuspend();
#endif
#endif
  CthResume(next);
}

void CthAwaken(CthThread th)
{
  if (B(th)->awakenfn == 0) CthNoStrategy();
#ifndef CMK_OPTIMIZE
#if ! CMK_TRACE_IN_CHARM
  if(CpvAccess(traceOn))
    traceAwaken(th);
#endif
#endif
  B(th)->awakenfn(B(th)->token, CQS_QUEUEING_FIFO, 0, 0);
	B(th)->scheduled = 1;
}

void CthYield()
{
  CthAwaken(CthCpvAccess(CthCurrent));
  CthSuspend();
}

void CthAwakenPrio(CthThread th, int s, int pb, unsigned int *prio)
{
  if (B(th)->awakenfn == 0) CthNoStrategy();
#ifndef CMK_OPTIMIZE
#if ! CMK_TRACE_IN_CHARM
  if(CpvAccess(traceOn))
    traceAwaken(th);
#endif
#endif
  B(th)->awakenfn(B(th)->token, s, pb, prio);
	B(th)->scheduled = 1;
}

void CthYieldPrio(int s, int pb, unsigned int *prio)
{
  CthAwakenPrio(CthCpvAccess(CthCurrent), s, pb, prio);
  CthSuspend();
}

/*
	Add a new thread listener to a thread 
*/
void CthAddListener(CthThread t,struct CthThreadListener *l){
		struct CthThreadListener *p=B(t)->listener;
		if(p== NULL){ /* first listener */
			B(t)->listener=l;
			l->thread = t;
			l->next=NULL;
			return;	
		}
		/* Add l at end of current chain of listeners: */
		while(p->next != NULL){
			p = p->next;
		}
		p->next = l;
		l->next = NULL;
		l->thread = t;
}

/*************************** Stack-Copying Threads (obsolete) *******************
Basic idea: switch from thread A (currently running) to thread B by copying
A's stack from the system stack area into A's buffer in the heap, then
copy B's stack from its heap buffer onto the system stack.

This allows thread migration, because the system stack is in the same
location on every processor; but the context-switching overhead (especially
for threads with deep stacks) is extremely high.

Written by Josh Yelon around May 1999
*/
#if CMK_THREADS_USE_STACKCOPY

#define SWITCHBUF_SIZE 16384

typedef struct CthProcInfo *CthProcInfo;

typedef struct CthThreadStruct
{
  CthThreadBase base;
  CthVoidFn  startfn;    /* function that thread will execute */
  void      *startarg;   /* argument that start function will be passed */
  qt_t      *savedstack; /* pointer to saved stack */
  int        savedsize;  /* length of saved stack (zero when running) */	
  int        stacklen;   /* length of the allocated savedstack >= savedsize */
  qt_t      *savedptr;   /* stack pointer */	
} CthThreadStruct;

CthThread CthPup(pup_er p, CthThread t)
{
    if (pup_isUnpacking(p))
      { t = (CthThread) malloc(sizeof(CthThreadStruct));_MEMCHECK(t);}
    pup_bytes(p, (void*) t, sizeof(CthThreadStruct)); 
    CthPupBase(p,&t->base,0);
    pup_int(p,&t->savedsize);
    if (pup_isUnpacking(p)) {
      t->savedstack = (qt_t*) malloc(t->savedsize);_MEMCHECK(t->savedstack);
    }
    pup_bytes(p, (void*) t->savedstack, t->savedsize);
    
    if (pup_isDeleting(p))
      {CthFree(t);t=0;}

    return t;
}

struct CthProcInfo
{
  qt_t      *stackbase;
  qt_t      *switchbuf_sp;
  qt_t      *switchbuf;
};

CthCpvDeclare(CthProcInfo, CthProc);

static void CthThreadInit(CthThread t, CthVoidFn fn, void *arg)
{
  CthThreadBaseInit(&t->base);
  t->startfn = fn;
  t->startarg = arg;
  t->savedstack = 0;
  t->savedsize = 0;
  t->stacklen = 0;
  t->savedptr = 0;
}

static void CthThreadFree(CthThread t)
{
  CthThreadBaseFree(&t->base);
  if (t->savedstack) free(t->savedstack);
  free(t);
}

void CthFree(t)
CthThread t;
{
  if (t==NULL) return;
  CthProcInfo proc = CthCpvAccess(CthProc);	

  if (t != CthSelf()) {
    CthThreadFree(t);
  } else
    t->base.exiting = 1;
}

void CthDummy() { }

void CthInit(char **argv)
{
  CthThread t; CthProcInfo p; qt_t *switchbuf, *sp;

  CthCpvInitialize(CthProcInfo, CthProc);

  CthBaseInit(argv);
  t = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(t);
  CthCpvAccess(CthCurrent)=t;
  CthThreadInit(t,0,0);

  p = (CthProcInfo)malloc(sizeof(struct CthProcInfo));
  _MEMCHECK(p);
  CthCpvAccess(CthProc)=p;
  /* leave some space for current stack frame < 256 bytes */
  sp = (qt_t*)(((size_t)&t) & ~((size_t)0xFF));
  p->stackbase = QT_SP(sp, 0x100);
  switchbuf = (qt_t*)malloc(QT_STKALIGN + SWITCHBUF_SIZE);
  _MEMCHECK(switchbuf);
  switchbuf = (qt_t*)((((size_t)switchbuf)+QT_STKALIGN) & ~(QT_STKALIGN-1));
  p->switchbuf = switchbuf;
  sp = QT_SP(switchbuf, SWITCHBUF_SIZE);
  sp = QT_ARGS(sp,0,0,0,(qt_only_t*)CthDummy);
  p->switchbuf_sp = sp;

  CmiThreadIs_flag |= CMI_THREAD_IS_STACKCOPY;
}

static void CthOnly(CthThread t, void *dum1, void *dum2)
{
  t->startfn(t->startarg);
  CthThreadFinished(t);
}

static void CthResume1(qt_t *sp, CthProcInfo proc, CthThread t)
{
  int bytes; qt_t *lo, *hi;
  CthThread old = CthCpvAccess(CthCurrent);
  CthBaseResume(t);
  if (old->base.exiting) {
    CthThreadFree(old);
  } else {
#ifdef QT_GROW_DOWN
    lo = sp; hi = proc->stackbase;
#else
    hi = sp; lo = proc->stackbase;
#endif
    bytes = ((size_t)hi)-((size_t)lo);
    if(bytes > old->stacklen) {
      if(old->savedstack) free((void *)old->savedstack);
      old->savedstack = (qt_t*)malloc(bytes);
      _MEMCHECK(old->savedstack);
      old->stacklen = bytes;
    }
    old->savedsize = bytes;
    old->savedptr = sp;
    memcpy(old->savedstack, lo, bytes);
  }
  if (t->savedstack) {
#ifdef QT_GROW_DOWN
    lo = t->savedptr;
#else
    lo = proc->stackbase;
#endif
    memcpy(lo, t->savedstack, t->savedsize);
    t->savedsize=0;
    sp = t->savedptr;
  } else {
    sp = proc->stackbase;
    sp = QT_ARGS(sp,t,0,0,(qt_only_t*)CthOnly);
  }
  QT_ABORT((qt_helper_t*)CthDummy,0,0,sp);
}

void CthResume(t)
CthThread t;
{
  CthProcInfo proc = CthCpvAccess(CthProc);
  QT_BLOCK((qt_helper_t*)CthResume1, proc, t, proc->switchbuf_sp);
}

CthThread CthCreate(CthVoidFn fn,void *arg,int size)
{
  CthThread result = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(result);
  CthThreadInit(result, fn, arg);
  return result;
}
CthThread CthCreateMigratable(CthVoidFn fn,void *arg,int size)
{
  /*All threads are migratable under stack copying*/
  return CthCreate(fn,arg,size);
}

/**************************************************************************
QuickThreads does not work on Win32-- our stack-shifting large allocas
fail a stack depth check.  Windows NT and 98 provide a user-level thread 
interface called "Fibers", used here.

Written by Sameer Paranjpye around October 2000
*/
#elif  CMK_THREADS_ARE_WIN32_FIBERS
#include <windows.h>
#include <winbase.h>

#ifndef _WIN32_WINNT
#define _WIN32_WINNT  0x0400
#endif

#if(_WIN32_WINNT >= 0x0400)
typedef VOID (WINAPI *PFIBER_START_ROUTINE)(
    LPVOID lpFiberParameter
    );
typedef PFIBER_START_ROUTINE LPFIBER_START_ROUTINE;
#endif

#if(_WIN32_WINNT >= 0x0400)
WINBASEAPI
LPVOID
WINAPI
CreateFiber(
    DWORD dwStackSize,
    LPFIBER_START_ROUTINE lpStartAddress,
    LPVOID lpParameter
    );

WINBASEAPI
VOID
WINAPI
DeleteFiber(
    LPVOID lpFiber
    );

WINBASEAPI
LPVOID
WINAPI
ConvertThreadToFiber(
    LPVOID lpParameter
    );

WINBASEAPI
VOID
WINAPI
SwitchToFiber(
    LPVOID lpFiber
    );

WINBASEAPI
BOOL
WINAPI
SwitchToThread(
    VOID
    );
#endif /* _WIN32_WINNT >= 0x0400 */


struct CthThreadStruct
{
  CthThreadBase base;
  LPVOID     fiber;
};

CthCpvStatic(CthThread,  CthPrevious);

typedef CthThread threadTable[100];
CthCpvStatic(threadTable, exitThreads);
CthCpvStatic(int, 	nExit);

static void CthThreadInit(CthThread t)
{
  CthThreadBaseInit(&t->base);
}

void CthInit(char **argv)
{
  CthThread t;

  CthCpvInitialize(CthThread,  CthPrevious);
  CthCpvInitialize(int,        nExit);
  CthCpvInitialize(threadTable,        exitThreads);

  CthCpvAccess(CthPrevious)=0;
  CthCpvAccess(nExit)=0;

  CthBaseInit(argv);  
  t = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(t);
  CthCpvAccess(CthCurrent)=t;
  CthThreadInit(t);
  t->fiber = ConvertThreadToFiber(t);
  _MEMCHECK(t->fiber);

  CmiThreadIs_flag |= CMI_THREAD_IS_FIBERS;
}

void CthThreadFree(CthThread old)
{
  CthThreadBaseFree(&old->base);
  if (old->fiber) DeleteFiber((PVOID)old->fiber);
  free(old);
}

static void CthClearThreads()
{
  int i,p,m;
  int n = CthCpvAccess(nExit);
  CthThread tc = CthCpvAccess(CthCurrent);
  CthThread tp = CthCpvAccess(CthPrevious);
  m = n;
  p=0;
  for (i=0; i<m; i++) {
    CthThread t = CthCpvAccess(exitThreads)[i];
    if (t && t != tc && t != tp) {
      CthThreadFree(t);
      CthCpvAccess(nExit) --;
    }
    else {
      if (p != i) CthCpvAccess(exitThreads)[p] = t;
      p++;
    }
  }
  if (m!=p)
  for (i=m; i<n; i++,p++) {
      CthCpvAccess(exitThreads)[p] = CthCpvAccess(exitThreads)[i];
  }
}

void CthFree(CthThread t)
{
  if (t==NULL) return;
  /* store into exiting threads table to avoid delete thread itself */
  CthCpvAccess(exitThreads)[CthCpvAccess(nExit)++] = t;
  if (t==CthCpvAccess(CthCurrent)) 
  {
     t->base.exiting = 1;
  } 
  else 
  {
    CthClearThreads();
/*  was
    if (t->data) free(t->data);
    DeleteFiber(t->fiber);
    free(t);
*/
  }
}

#if 0
void CthFiberBlock(CthThread t)
{
  CthThread tp;
  
  SwitchToFiber(t->fiber);
  tp = CthCpvAccess(CthPrevious);
  if (tp != 0 && tp->killed == 1)
    CthThreadFree(tp);
}
#endif

void CthResume(CthThread t)
{
  CthThread tc;

  tc = CthCpvAccess(CthCurrent);
  if (t == tc) return;
  CthBaseResume(t);
  CthCpvAccess(CthPrevious)=tc;
#if 0
  if (tc->base.exiting) 
  {
    SwitchToFiber(t->fiber);
  } 
  else 
    CthFiberBlock(t);
#endif
  SwitchToFiber(t->fiber);
}

VOID CALLBACK FiberSetUp(PVOID fiberData)
{
  void **ptr = (void **) fiberData;
  qt_userf_t* fn = (qt_userf_t *)ptr[0];
  void *arg = ptr[1];
  CthThread  t = CthSelf();

  CthClearThreads();

  fn(arg);

  CthCpvAccess(exitThreads)[CthCpvAccess(nExit)++] = t;
  CthThreadFinished(t);
}

CthThread CthCreate(CthVoidFn fn, void *arg, int size)
{
  CthThread result; 
  void**    fiberData;
  fiberData = (void**)malloc(2*sizeof(void *));
  fiberData[0] = (void *)fn;
  fiberData[1] = arg;
  
  result = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(result);
  CthThreadInit(result);
  result->fiber = CreateFiber(size, FiberSetUp, (PVOID) fiberData);
  if (!result->fiber)
    CmiAbort("CthCreate failed to create fiber!\n");
  
  return result;
}

CthThread CthPup(pup_er p, CthThread t)
{
  CmiAbort("CthPup not implemented.\n");
  return 0;
}
CthThread CthCreateMigratable(CthVoidFn fn,void *arg,int size)
{
  /*Fibers are never migratable, unless we can figure out how to set their stacks*/
  return CthCreate(fn,arg,size);
}

/***************************************************
Use Posix Threads to simulate cooperative user-level
threads.  This version is very portable but inefficient.

Written by Milind Bhandarkar around November 2000
*/
#elif CMK_THREADS_USE_PTHREADS

#include <pthread.h>
#include <errno.h> /* for perror */

struct CthThreadStruct
{
  CthThreadBase base;
  pthread_t  self;
  pthread_cond_t cond;
  pthread_cond_t *creator;
  CthVoidFn  fn;
  void      *arg;
  char       inited;
};

/**
  The sched_mutex is the current token of execution.
  Only the running thread holds this lock; all other threads
  have released the lock and are waiting on their condition variable.
*/
CthCpvStatic(pthread_mutex_t, sched_mutex);

static void CthThreadInit(t)
CthThread t;
{
  CthThreadBaseInit(&t->base);
  t->inited = 0;
  pthread_cond_init(&(t->cond) , (pthread_condattr_t *) 0);
}

void CthInit(char **argv)
{
  CthThread t;

  CthCpvInitialize(pthread_mutex_t, sched_mutex);

  pthread_mutex_init(&CthCpvAccess(sched_mutex), (pthread_mutexattr_t *) 0);
  pthread_mutex_lock(&CthCpvAccess(sched_mutex));
  CthBaseInit(argv); 
  t = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(t);
  CthCpvAccess(CthCurrent)=t;
  CthThreadInit(t);
  t->self = pthread_self();

  CmiThreadIs_flag |= CMI_THREAD_IS_PTHREADS;
}

void CthFree(t)
CthThread t;
{
  if (t==NULL) return;
  if (t==CthCpvAccess(CthCurrent)) {
    t->base.exiting = 1;
  } else {
    CthThreadBaseFree(&t->base);
    free(t);
  }
}

void CthResume(CthThread t)
{
  CthThread tc = CthCpvAccess(CthCurrent);
  if (t == tc) return;
  CthBaseResume(t);
  pthread_cond_signal(&(t->cond)); /* wake up the next thread */
  if (tc->base.exiting) {
    pthread_mutex_unlock(&CthCpvAccess(sched_mutex));
    pthread_exit(0);
  } else {
    /* pthread_cond_wait might (with low probability) return when the 
      condition variable has not been signaled, guarded with 
      predicate checks */
    do {
    pthread_cond_wait(&(tc->cond), &CthCpvAccess(sched_mutex));
    } while (tc!=CthCpvAccess(CthCurrent)) ;
  }
}

static void *CthOnly(void * arg)
{
  CthThread th = (CthThread)arg;
  th->inited = 1;
  pthread_detach(pthread_self());
  pthread_mutex_lock(&CthCpvAccess(sched_mutex));
  pthread_cond_signal(th->creator);
  do {
  pthread_cond_wait(&(th->cond), &CthCpvAccess(sched_mutex));
  } while (arg!=CthCpvAccess(CthCurrent)) ;
  th->fn(th->arg);
  CthThreadFinished(th);
  return 0;
}

CthThread CthCreate(CthVoidFn fn, void *arg, int size)
{
  static int reported = 0;
  pthread_attr_t attr;
  int r;
  CthThread result;
  CthThread self = CthSelf();
  /* size is ignored in this version */
  result = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(result);
  CthThreadInit(result);
  result->fn = fn;
  result->arg = arg;
  result->creator = &(self->cond);

  /* try set pthread stack, not necessarily supported on all platforms */
  pthread_attr_init(&attr);
  if (size<1024) size = CthCpvAccess(_defaultStackSize);
  if (0!=(r=pthread_attr_setstacksize(&attr,size))) {
      if (!reported) {
        CmiPrintf("Warning: pthread_attr_setstacksize failed");
	errno = r;
	perror("pthread_attr_setstacksize");
        reported = 1;
      }
  }

  if (0 != pthread_create(&(result->self), &attr, CthOnly, (void*) result)) 
    CmiAbort("CthCreate failed to created a new pthread\n");
  do {
  pthread_cond_wait(&(self->cond), &CthCpvAccess(sched_mutex));
  } while (result->inited==0);
  return result;
}

CthThread CthPup(pup_er p, CthThread t)
{
  CmiAbort("CthPup not implemented.\n");
  return 0;
}
CthThread CthCreateMigratable(CthVoidFn fn,void *arg,int size)
{
  /*Pthreads are never migratable, unless we can figure out how to set their stacks*/
  return CthCreate(fn,arg,size);
}

/***************************************************************
Use SysV r3 setcontext/getcontext calls instead of
quickthreads.  This works on lots of architectures (such as
SUN, IBM SP, O2K, DEC Alpha, IA64, Cray X1, Linux with newer version of 
 glibc such as the one with RH9) 
On some machine such as IA64 and Cray X1, the context version is the 
only thread package that is working. 

Porting should be easy. To port context threads, one need to set the 
direction of the thread stack properly in conv-mach.h.

Note: on some machine like Sun and IBM SP, one needs to link with memory gnuold
 to have this context thread working.

Written by Gengbin Zheng around April 2001
*/
#elif (CMK_THREADS_USE_CONTEXT || CMK_THREADS_USE_JCONTEXT)

#include <signal.h>
#include <errno.h>

#if CMK_THREADS_USE_CONTEXT
/* system builtin context routines: */
#include <ucontext.h>

#define uJcontext_t ucontext_t
#define setJcontext setcontext
#define getJcontext getcontext
#define swapJcontext swapcontext
#define makeJcontext makecontext
typedef void (*uJcontext_fn_t)(void);

#else /* CMK_THREADS_USE_JCONTEXT */
/* Orion's setjmp-based context routines: */
#include <uJcontext.h>
#include <uJcontext.c>

#endif


struct CthThreadStruct
{
  CthThreadBase base;
  double * dummy;
  uJcontext_t context;
};


static void CthThreadInit(t)
CthThread t;
{
  CthThreadBaseInit(&t->base);
}

/* Threads waiting to be destroyed */
CpvStaticDeclare(CthThread , doomedThreadPool);

void CthInit(char **argv)
{
  CthThread t;

  CthBaseInit(argv);
  t = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(t);
  CthCpvAccess(CthCurrent)=t;
  if (0 != getJcontext(&t->context))
    CmiAbort("CthInit: getcontext failed.\n");
  CthThreadInit(t);
  CpvInitialize(CthThread, doomedThreadPool);
  CpvAccess(doomedThreadPool) = (CthThread)NULL;

  /* don't trust the _defaultStackSize */
#ifdef MINSIGSTKSZ
  if (CthCpvAccess(_defaultStackSize) < MINSIGSTKSZ) 
    CthCpvAccess(_defaultStackSize) = MINSIGSTKSZ;
#endif
#if CMK_THREADS_USE_CONTEXT
  CmiThreadIs_flag |= CMI_THREAD_IS_CONTEXT;
#else
  CmiThreadIs_flag |= CMI_THREAD_IS_UJCONTEXT;
#endif
}

static void CthThreadFree(CthThread t)
{
  /* avoid freeing thread while it is being used, store in pool and 
     free it next time. Note the last thread in pool won't be free'd! */
  CthThread doomed=CpvAccess(doomedThreadPool);
  CpvAccess(doomedThreadPool) = t;
  if (doomed != NULL) {
    CthThreadBaseFree(&doomed->base);
    free(doomed);
  }
}

void CthFree(CthThread t)
{
  if (t==NULL) return;
  if (t==CthCpvAccess(CthCurrent)) {
    t->base.exiting = 1; /* thread is already running-- free on next swap out */
  } else {
    CthThreadFree(t);
  }
}

void CthResume(t)
CthThread t;
{
  CthThread tc;
  tc = CthCpvAccess(CthCurrent);
  if (t != tc) { /* Actually switch threads */
    CthBaseResume(t);
    if (!tc->base.exiting) 
    {
      if (0 != swapJcontext(&tc->context, &t->context)) 
        CmiAbort("CthResume: swapcontext failed.\n");
    } 
    else /* tc->base.exiting, so jump directly to next context */ 
    {
      CthThreadFree(tc);
      setJcontext(&t->context);
    }
  }
/*This check will mistakenly fail if the thread migrates (changing tc)
  if (tc!=CthCpvAccess(CthCurrent)) { CmiAbort("Stack corrupted?\n"); }
*/
}

void CthStartThread(qt_userf_t fn,void *arg)
{
  fn(arg);
  CthThreadFinished(CthSelf());
}

#define STP_STKALIGN(sp, alignment) \
  ((void *)((((qt_word_t)(sp)) + (alignment) - 1) & ~((alignment)-1)))

int ptrDiffLen(const void *a,const void *b) {
	char *ac=(char *)a, *bc=(char *)b;
	int ret=ac-bc;
	if (ret<0) ret=-ret;
	return ret;
}

static CthThread CthCreateInner(CthVoidFn fn,void *arg,int size,int migratable)
{
  CthThread result;
  char *stack, *ss_sp, *ss_end;
  result = (CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(result);
  CthThreadInit(result);
#ifdef MINSIGSTKSZ
  if (size<MINSIGSTKSZ) size = CthCpvAccess(_defaultStackSize);
#endif
  CthAllocateStack(&result->base,&size,migratable);
  stack = result->base.stack;
  
  if (0 != getJcontext(&result->context))
    CmiAbort("CthCreateInner: getcontext failed.\n");

  ss_end = stack + size;

/**
 Decide where to point the uc_stack.ss_sp field of our "context"
 structure.  The configuration values CMK_CONTEXT_STACKBEGIN, 
 CMK_CONTEXT_STACKEND, and CMK_CONTEXT_STACKMIDDLE determine where to
 point ss_sp: to the beginning, end, and middle of the stack buffer
 respectively.  The default, used by most machines, is CMK_CONTEXT_STACKBEGIN.
*/
#if CMK_THREADS_USE_JCONTEXT /* Jcontext is always STACKBEGIN */
  ss_sp = stack;
#elif CMK_CONTEXT_STACKEND /* ss_sp should point to *end* of buffer */
  ss_sp = stack+size-MINSIGSTKSZ; /* the MINSIGSTKSZ seems like a hack */
  ss_end = stack;
#elif CMK_CONTEXT_STACKMIDDLE /* ss_sp should point to *middle* of buffer */
  ss_sp = stack+size/2;
#else /* CMK_CONTEXT_STACKBEGIN, the usual case  */
  ss_sp = stack;
#endif
  
  result->context.uc_stack.ss_sp = STP_STKALIGN(ss_sp,sizeof(char *)*8);
  result->context.uc_stack.ss_size = ptrDiffLen(result->context.uc_stack.ss_sp,ss_end);
  result->context.uc_stack.ss_flags = 0;
  result->context.uc_link = 0;
  
  CthAliasEnable(B(result)); /* Change to new thread's stack while building context */
  errno = 0;
  makeJcontext(&result->context, (uJcontext_fn_t)CthStartThread, 2, (void *)fn,(void *)arg);
  if(errno !=0) { 
    perror("makecontext"); 
    CmiAbort("CthCreateInner: makecontext failed.\n");
  }
  CthAliasEnable(B(CthCpvAccess(CthCurrent)));
  return result;  
}

CthThread CthCreate(CthVoidFn fn,void *arg,int size)
{
  return CthCreateInner(fn,arg,size,0);
}
CthThread CthCreateMigratable(CthVoidFn fn,void *arg,int size)
{
  return CthCreateInner(fn,arg,size,1);
}

CthThread CthPup(pup_er p, CthThread t)
{
  if (pup_isUnpacking(p)) {
	  t=(CthThread)malloc(sizeof(struct CthThreadStruct));
	  _MEMCHECK(t);
	  CthThreadInit(t);
  }
  CthPupBase(p,&t->base,1);
  
  /*Pup the processor context as bytes-- this is not guarenteed to work!*/
  pup_bytes(p,&t->context,sizeof(t->context));
  if (pup_isDeleting(p)) {
	  CthFree(t);
	  return 0;
  }
  return t;
}

#else 
/***************************************************************
Basic qthreads implementation. 

These threads can also add a "protection block" of
inaccessible memory to detect stack overflows, which
would otherwise just trash the heap.

(7/13/2001 creation times on 300MHz AMD K6-3 x86, Linux 2.2.18:
Qt setjmp, without stackprotect: 18.5 us
Qt i386, without stackprotect: 17.9 us
Qt setjmp, with stackprotect: 68.6 us
)

Written by Josh Yelon around 1995
*/

#if defined(CMK_OPTIMIZE) || (!CMK_MEMORY_PROTECTABLE)
#  define CMK_STACKPROTECT 0

#  define CthMemAlign(x,n) 0
#  define CthMemoryProtect(m,p,l) CmiAbort("Shouldn't call CthMemoryProtect!\n")
#  define CthMemoryUnprotect(m,p,l) CmiAbort("Shouldn't call CthMemoryUnprotect!\n")
#else
#  define CMK_STACKPROTECT 1

  extern void setProtection(char*, char*, int, int);
#  include "sys/mman.h"
#  define CthMemAlign(x,n) memalign((x),(n))
#  define CthMemoryProtect(m,p,l) mprotect(p,l,PROT_NONE);setProtection((char*)m,p,l,1);
#  define CthMemoryUnprotect(m,p,l) mprotect(p,l,PROT_READ | PROT_WRITE);setProtection((char*)m,p,l,0);
#endif

struct CthThreadStruct
{
  CthThreadBase base;

  char      *protect;
  int        protlen;

  qt_t      *stack;
  qt_t      *stackp;
};

static CthThread CthThreadInit(void)
{
  CthThread ret=(CthThread)malloc(sizeof(struct CthThreadStruct));
  _MEMCHECK(ret);
  CthThreadBaseInit(&ret->base);
  ret->protect = 0;
  ret->protlen = 0;
  
  return ret;
}

static void CthThreadFree(CthThread t)
{
  CthThreadBaseFree(&t->base);
  if (t->protlen!=0) {
    CthMemoryUnprotect(t->stack, t->protect, t->protlen);
    free(t->stack);
  }
  free(t);
}

void CthInit(char **argv)
{
  CthThread mainThread;

  CthBaseInit(argv);  
  mainThread=CthThreadInit();
  CthCpvAccess(CthCurrent)=mainThread;
  /* mainThread->base.suspendable=0;*/ /*Can't suspend main thread (trashes Quickthreads jump buffer)*/

#if CMK_THREADS_ALIAS_STACK
  CmiThreadIs_flag |= CMI_THREAD_IS_QT_ALIAS;
#else
  CmiThreadIs_flag |= CMI_THREAD_IS_QT;
#endif
}

void CthFree(CthThread t)
{
  if (t==NULL) return;
  if (t==CthCpvAccess(CthCurrent)) {
    t->base.exiting = 1;
  } else {
    CthThreadFree(t);
  }
}

static void *CthAbortHelp(qt_t *sp, CthThread old, void *null)
{
  CthThreadFree(old);
  return (void *) 0;
}

static void *CthBlockHelp(qt_t *sp, CthThread old, void *null)
{
  old->stackp = sp;
  return (void *) 0;
}

void CthResume(t)
CthThread t;
{
  CthThread tc = CthCpvAccess(CthCurrent);
  if (t == tc) return;
  CthBaseResume(t);
  if (tc->base.exiting) {
    QT_ABORT((qt_helper_t*)CthAbortHelp, tc, 0, t->stackp);
  } else {
    QT_BLOCK((qt_helper_t*)CthBlockHelp, tc, 0, t->stackp);
  }
/*This check will mistakenly fail if the thread migrates (changing tc)
  if (tc!=CthCpvAccess(CthCurrent)) { CmiAbort("Stack corrupted?\n"); }
*/
}

static void CthOnly(void *arg, void *vt, qt_userf_t fn)
{
  fn(arg);
  CthThreadFinished(CthSelf());
}

static CthThread CthCreateInner(CthVoidFn fn, void *arg, int size,int Migratable)
{
  CthThread result; qt_t *stack, *stackbase, *stackp;
  int doProtect=(!Migratable) && CMK_STACKPROTECT;
  result=CthThreadInit();
  if (doProtect) 
  { /*Can only protect on a page boundary-- allocate an extra page and align stack*/
	  if (size==0) size=CthCpvAccess(_defaultStackSize);
	  size = (size+(CMK_MEMORY_PAGESIZE*2)-1) & ~(CMK_MEMORY_PAGESIZE-1);
	  stack = (qt_t*)CthMemAlign(CMK_MEMORY_PAGESIZE, size);
  } else
	  stack=CthAllocateStack(&result->base,&size,Migratable);
  CthAliasEnable(B(result)); /* Change to new thread's stack while setting args */
  stackbase = QT_SP(stack, size);
  stackp = QT_ARGS(stackbase, arg, result, (qt_userf_t *)fn, CthOnly);
  CthAliasEnable(B(CthCpvAccess(CthCurrent)));
  result->stack = stack;
  result->stackp = stackp;
  if (doProtect) {
#ifdef QT_GROW_UP
  /*Stack grows up-- protect highest page of stack*/
    result->protect = ((char*)stack) + size - CMK_MEMORY_PAGESIZE;
#else
  /*Stack grows down-- protect lowest page in stack*/
    result->protect = ((char*)stack);
#endif
    result->protlen = CMK_MEMORY_PAGESIZE;
    CthMemoryProtect(stack, result->protect, result->protlen);
  }
  return result;
}
CthThread CthCreate(CthVoidFn fn, void *arg, int size)
{ return CthCreateInner(fn,arg,size,0);}

CthThread CthCreateMigratable(CthVoidFn fn, void *arg, int size)
{ return CthCreateInner(fn,arg,size,1);}

CthThread CthPup(pup_er p, CthThread t)
{
  if (pup_isUnpacking(p)) {
	  t=CthThreadInit();
  }
#if CMK_THREADS_ALIAS_STACK
  CthPupBase(p,&t->base,0);
#else
  CthPupBase(p,&t->base,1);
#endif
  
  /*Pup the stack pointer as bytes-- this works because stack is migratable*/
  pup_bytes(p,&t->stackp,sizeof(t->stackp));

  /*Don't worry about stack protection on migration*/  

  if (pup_isDeleting(p)) {
	  CthFree(t);
	  return 0;
  }
  return t;
}

#endif

void CthTraceResume(CthThread t)
{
  traceResume(&t->base.tid);
}

