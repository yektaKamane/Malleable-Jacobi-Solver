/// Base class to represent user object
/** This is what the user class becomes. It adds minimal functionality
    to the user's code, mostly OVT, and links back to the parent sim object. 
    It also provides the derived templated class chpt for checkpointing. */
#ifndef REP_H
#define REP_H
#include "sim.decl.h"

extern CProxy_sim POSE_Objects_RO;

// our configuration bundle
extern POSE_Config pose_config;

/// Base representation class
class rep 
{
 protected:
  /// Pointer to synchronization strategy
  strat *myStrat;          
 public:
#if POSE_COMM_ON
  CProxy_sim POSE_Objects;
#endif
  /// Pointer to poser wrapper
  sim *parent;             
  /// the object's unique handle
  /** Initialized to index of poser wrapper in POSE_objects array */
  int myHandle;            
  /// The object's virtual time (OVT)
  POSE_TimeType ovt;
  /// The object's real time (ORT)
  double ort;
  /// Checkpointed seed for POSE_rand
  unsigned int prand_seed;
  /// Checkpointed 48-bit seed to support uniform and linear rand
  unsigned short int prand48_seed[3];
  /// Flag indicating this object uses anti-methods rather than checkpoints
  int anti_methods;
#ifdef MEM_TEMPORAL
  TimePool *localTimePool;
#endif
  /// Basic Constructor
  rep() { 
    ovt = 0; ort = 0.0; parent = NULL; myStrat = NULL; 
    anti_methods = 0;
  #ifndef SEQUENTIAL_POSE
  #ifdef POSE_COMM_ON    
    POSE_Objects = POSE_Objects_RO;
    ComlibDelegateProxy(&POSE_Objects);
  #endif
  #endif
#ifdef MEM_TEMPORAL    
    localTimePool = (TimePool *)CkLocalBranch(TempMemID);
#endif
  }
  /// Initializing Constructor
  rep(POSE_TimeType init_ovt) { 
    ovt = init_ovt; ort = 0.0; anti_methods = 0;
      }
  /// Destructor
  virtual ~rep() {
     }
  /// Initializer called from poser wrapper constructor
  void init(eventMsg *m);
  /// Return the OVT
  POSE_TimeType OVT() { return ovt; }
  /// Set the OVT to t
  void SetOVT(POSE_TimeType t) { ovt = t; }
  /// Make object use anti-methods rather than checkpointing
  void useAntimethods() { anti_methods = 1; }
  /// Make object use checkpointing rather than anti-methods
  void turnOffAntimethods() { anti_methods = 0; }
  /// Check if this object uses anti-methods rather than checkpointing
  int usesAntimethods() { return anti_methods; }
  /// Elapse time by incrementing the OVT by dt
  void elapse(POSE_TimeType dt) { ovt += dt; }
  /// Update the OVT and ORT at event start to auto-elapse to event timestamp
  /** If event has timestamp > OVT, OVT elapses to timestamp, otherwise
      there is no change to OVT. ORT updates similarly. */
  void update(POSE_TimeType t, double rt);
  /// Called on every object at end of simulation
  virtual void terminus() { 
    //CkPrintf("Object %d terminus at time %d\n", myHandle, ovt);
  }
  /// Timestamps event message, sets priority, and records in spawned list
  virtual void registerTimestamp(int idx, eventMsg *m, POSE_TimeType offset);
  /// Assignment operator
  /** Derived classes must provide assignment */
  virtual rep& operator=(const rep& obj) { 
    ovt = obj.ovt; 
    ort = obj.ort;
    myHandle = obj.myHandle;
    anti_methods = obj.anti_methods;
    prand_seed = obj.prand_seed;
    prand48_seed[0]=obj.prand48_seed[0];
    prand48_seed[1]=obj.prand48_seed[1];
    prand48_seed[2]=obj.prand48_seed[2];
    return *this;
  }
  /// Dump all data fields
#if USE_LONG_TIMESTAMPS
  virtual void dump() { CkPrintf("[REP: ovt=%lld]\n", ovt); }
#else
  virtual void dump() { CkPrintf("[REP: ovt=%d]\n", ovt); }
#endif
  /// Pack/unpack/sizing operator
  /** Derived classes must provide pup */
  virtual void pup(PUP::er &p) { 
    p(ovt); p(ort); p(myHandle); p(anti_methods); p(prand48_seed,3);
  }
#ifdef SEQUENTIAL_POSE
  void checkpoint(rep *) { }
  void restore(rep *) { }
#endif  
  void POSE_srand(unsigned int pseed) { prand_seed = pseed; }
  int POSE_rand() { 
    int rnum; 
    srand(prand_seed); 
    rnum = rand(); 
    prand_seed = rnum+INT_MAX; 
    return rnum; 
  }
  unsigned int POSE_urand() { 
    int rnum; 
    srand(prand_seed); 
    rnum = rand(); 
    prand_seed = rnum+INT_MAX; 
    return prand_seed;
  }
  inline long int POSE_Linear_rand() { return nrand48(prand48_seed); }
  inline double POSE_Uniform_rand() { return erand48(prand48_seed); }

};

#endif
