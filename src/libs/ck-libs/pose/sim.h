/// Sim is the base class for all poser entities
/** This is the wrapper class that encapsulates synchronization strategy,
    representation object and event queue, and manages the orchestration
    of these components. 
    This file also defines the three basic POSE messages: 
    eventMsg, from which all user messages inherit; 
    cancelMsg, for cancelling events; 
    and prioMsg, a null message.
    All three have a priority field which is set to the timestamp; thus all 
    messages with earlier timestamp have higher priority. */
#ifndef SIM_H
#define SIM_H
#include "sim.decl.h"
#include <stdarg.h>

extern CProxy_sim POSE_Objects; 
extern CkChareID POSE_Coordinator_ID; 
class sim; // needed for eventMsg definition below

/// All user event messages inherit from this
/** Adds timestamp and event ID to message, plus other info useful for the
    underlying simulation layer.  Prioritized by default, and given a priority
    based on the timestamp. Events which take no parameters must
    still pass an eventMsg. */
class eventMsg : public CMessage_eventMsg {
public:
  /// The event's timestamp
  int timestamp;  
  /// The event's globally unique ID
  eventID evID;   
  /// The message size, used for message recycling (currently not used)
  int msgSize;    
  /// The sending PE, used for comm-based LB
  int fromPE;     
  /// Pointer to a poser wrapper; used to send the pointer to rep object
  sim *parent;    
  /// Pointer to synchronization strategy; used when creating rep object
  strat *str;
  /// Relative start time: for computing degree of parallelization
  double rst;
  /// Basic Constructor
  eventMsg() { }
  /// Destructor
  virtual ~eventMsg() { }
  /// Timestamps this message and generates a unique event ID
  /** Timestamps this message and generates a unique event ID for the event
      to be invoked on the receiving side.  Sets the priority of this
      message to timestamp - INT_MAX. */
  void Timestamp(int t) { 
    timestamp = t;  evID = GetEventID();  
#ifdef PRIO_MSGS
    setPriority(t-INT_MAX); 
#endif
    rst = 0.0;
  }
  /// Assignment operator: copies priority too
  eventMsg& operator=(const eventMsg& obj) {
    timestamp = obj.timestamp;
    evID = obj.evID;
    parent = obj.parent;
    str = obj.str;
    msgSize = obj.msgSize;
    rst = obj.rst;
#ifdef PRIO_MSGS    
    setPriority(timestamp-INT_MAX);
#endif
    return *this;
  }
  /// Allocates event message with space for priority
  /** This can also handle event message recycling (currently off) */
  void *operator new (size_t size) {  
#ifdef MSG_RECYCLING
    MemoryPool *localPool = (MemoryPool *)CkLocalBranch(MemPoolID);
    if (localPool->CheckPool(size) > 0)
      return localPool->GetBlock(size);
    else {
#endif
#ifdef PRIO_MSGS
      void *msg = CkAllocMsg(CMessage_eventMsg::__idx, size, 8*sizeof(int));
#else
      void *msg = CkAllocMsg(CMessage_eventMsg::__idx, size, 0);
#endif
      ((eventMsg *)msg)->msgSize = size;
      return msg;
#ifdef MSG_RECYCLING
    }
#endif
  }
  void operator delete(void *p) { 
#ifdef MSG_RECYCLING
    MemoryPool *localPool = (MemoryPool *)CkLocalBranch(MemPoolID);
    int ps = localPool->CheckPool(((eventMsg *)p)->msgSize);
    if (0) { // ((ps < MAX_POOL_SIZE) && (ps > -1)) {
      size_t msgSize = ((eventMsg *)p)->msgSize;
      memset(p, 0, msgSize);
      ((eventMsg *)p)->msgSize = msgSize;
      localPool->PutBlock(msgSize, p);
    }
    else
#endif
      CkFreeMsg(p);
  }
  /// Set priority field and queuing strategy
  void setPriority(int prio) {  
    *((int*)CkPriorityPtr(this)) = prio;
    CkSetQueueing(this, CK_QUEUEING_IFIFO);
  }
};

/// Cancellation message
class cancelMsg : public CMessage_cancelMsg {
public:
  /// Event to cancel
  /** Only this is needed to find the event to cancel */
  eventID evID;
  /// Timestamp of event to be cancelled
  /** Providing this makes finding the event faster */
  int timestamp;          
  /// Allocate cancellation message with priority field
  void *operator new (size_t size) {  
    return CkAllocMsg(CMessage_cancelMsg::__idx, size, 8*sizeof(int));
  } 
  /// Delete cancellation message
  void operator delete(void *p) {  CkFreeMsg(p);  }
  /// Set priority field and queuing strategy
  void setPriority(int prio) {
    *((int*)CkPriorityPtr(this)) = prio;
    CkSetQueueing(this, CK_QUEUEING_IFIFO);
  }
};

/// Prioritized null msg; used to sort Step calls
class prioMsg : public CMessage_prioMsg {
public:
  /// Allocate prioritized message with priority field
  void *operator new (size_t size) {
    return CkAllocMsg(CMessage_eventMsg::__idx, size, 8*sizeof(int));
  }
  /// Delete prioritized message
  void operator delete(void *p) {  CkFreeMsg(p);  }
  /// Set priority field and queuing strategy
  void setPriority(int prio) {
    *((int*)CkPriorityPtr(this)) = prio;
    CkSetQueueing(this, CK_QUEUEING_IFIFO);
  }
};

/// Used to specify a destination processor to migrate to during load balancing
class destMsg : public CMessage_destMsg {
public:
  int destPE;
};


/// Poser wrapper base class
/** The poser base class: all user posers are translated to classes that
    inherit from this class, and act as wrappers around the actual user 
    object's representation to control the simulation behavior.  These 
    objects are plugged into the POSE_objects array which is of this type. */
class sim : public ArrayElement1D {
 protected:
  /// Flag to indicate that a Step message is scheduled
  /** Need to re-evaluate the need/function of this... also, how is it used
      during load balancing... */
  int active; 
 public:
  /// This poser's event queue
  eventQueue *eq;
  /// This poser's synchronization strategy   
  strat *myStrat;
  /// This poser's user representation
  rep *objID;       
  /// List of incoming cancellations for this poser
  CancelList cancels;
  /// The local PVT to report to
  PVT *localPVT;    
  /// Unique global ID for this object on PVT branch
  int myPVTidx;
  /// Unique global ID for this object in load balancing data structures
  int myLBidx;
  /// Number of forward execution steps
  /** Is this needed/used? This is load balancing data... */
  int DOs;
  /// Number of undone events
  /** Is this needed/used? This is load balancing data... */
  int UNDOs;
  /// Synchronization strategy type (optimistic or conservative)
  int sync;
  /// Number of sends/recvs per PE
  int *srVector;    
  /// Most recent GVT estimate
  int lastGVT;
  /// Relative start time, start time, end time and current time
  /** Used to calculate degree of parallelism */
  double st, et, ct;
#ifdef POSE_STATS_ON
  /// The local statistics collector
  localStat *localStats; 
#endif
#ifdef LB_ON
  /// The local load balancer
  LBgroup *localLBG;
#endif
  /// Basic Constructor
  sim(void);
  sim(CkMigrateMessage *) {};
  /// Destructor
  virtual ~sim();
  /// Pack/unpack/sizing operator
  virtual void pup(PUP::er &p) {
    ArrayElement1D::pup(p); // call parent class pup method
    // pup simple types
    p(active); p(myPVTidx); p(myLBidx); p(sync); p(DOs); p(UNDOs);
    // pup event queue
    if (p.isUnpacking())
      eq = new eventQueue();
    eq->pup(p);
    // pup cancellations
    cancels.pup(p);
    if (p.isUnpacking()) { // reactivate migrated object
#ifdef POSE_STATS_ON
      localStats = (localStat *)CkLocalBranch(theLocalStats);
#endif
      localPVT = (PVT *)CkLocalBranch(ThePVT);
#ifdef LB_ON
      localLBG = TheLBG.ckLocalBranch();
#endif
      active = 0;
      myPVTidx = localPVT->objRegister(thisIndex, localPVT->getGVT(), sync, this);
#ifdef LB_ON
      myLBidx = localLBG->objRegister(thisIndex, sync, this);
#endif
    }
    else if (p.isPacking()) { // deactivate migrating object
      active = -1;
      localPVT->objRemove(myPVTidx);
#ifdef LB_ON
      localLBG->objRemove(myLBidx);
#endif
    }
  }
  /// Start a forward execution step on myStrat
  void Step();                 
  /// Start a prioritized forward execution step on myStrat
  void Step(prioMsg *m);       
  /// Report safe time to PVT branch
  void Status() { localPVT->objUpdate(myPVTidx, myStrat->SafeTime(), -1, -1); }
  /// Commit events based on new GVT estimate
  void Commit();               
  /// Add m to cancellation list
  void Cancel(cancelMsg *m); 
  /// Report load information to local load balancer
  void ReportLBdata();
  /// Migrate this poser to processor indicated in m
  void Migrate(destMsg *m) { migrateMe(m->destPE); }
  /// Return this poser's unique index on PVT branch
  int PVTindex() { return myPVTidx; }
  /// Test active flag
  int IsActive() { return active; }
  /// Set active flag
  void Activate() { active = 1; }
  /// Unset active flag
  void Deactivate() { active = 0; }
  /// Invoke an event on this poser according to fnIdx and pass it msg
  /** ResolveFn is generated along with the rest of the wrapper object and
      should handle all possible events on a poser. */
  virtual void ResolveFn(int fnIdx, void *msg) { }
  /// Invoke the commit version of an event to handle special behaviors
  /** This invokes the <fn>_commit method that user provides.  It can be
      used to perform special activities, statistics gathering, output, or
      whatever the user wishes. */
  virtual void ResolveCommitFn(int fnIdx, void *msg) { }
  /// Notify the PVT of a message send
  void registerSent(int timestamp) {
    localPVT->objUpdate(timestamp, SEND);
  }
  /// Used for buffered output
  /** Output is only printed when the event is committed */
  void CommitPrintf(const char *Fmt, ...) {
    va_list ap;
    va_start(ap,Fmt);
    InternalCommitPrintf(Fmt, ap);
    va_end(ap);
  }
  /// Used for buffered output of error messages
  /** Output is only printed when the event is committed */
  void CommitError(const char *Fmt, ...) {
    va_list ap;
    va_start(ap,Fmt);
    InternalCommitPrintf(Fmt, ap);
    va_end(ap);
    myStrat->currentEvent->commitErr = 1;
  }
  /// Dump all data fields
  void dump();
 private:
  /// Used by buffered print functions
  void InternalCommitPrintf (const char *Fmt, va_list ap) {
    char *tmp;
    size_t tmplen=myStrat->currentEvent->commitBfrLen + strlen(Fmt) + 1 +100;
    if (!(tmp = (char *)malloc(tmplen * sizeof(char)))) {
      CkPrintf("ERROR: sim::CommitPrintf: OUT OF MEMORY!\n");
      CkExit();
    }
    if (myStrat->currentEvent->commitBfr && myStrat->currentEvent->commitBfrLen) {
      strcpy(tmp, myStrat->currentEvent->commitBfr);
      free(myStrat->currentEvent->commitBfr);
      vsnprintf(tmp+strlen(tmp), tmplen, Fmt, ap); 
    }
    else vsnprintf(tmp, tmplen, Fmt, ap); 
    myStrat->currentEvent->commitBfrLen = strlen(tmp) + 1;  
    myStrat->currentEvent->commitBfr = tmp;
  }
};

#endif









