/// Global Virtual Time estimation for POSE
/** Implements the Global Virtual Time (GVT) algorithm; provides chare
    groups PVT and GVT. Objects interact with the local PVT branch.
    PVT branches summarize object info and report to the single
    floating GVT object, which broadcasts results to all PVT branches. */

#ifndef GVT_H
#define GVT_H
#include "gvt.decl.h"

/// synchronization strategies
#define OPTIMISTIC 0
#define CONSERVATIVE 1

/// Global handles ThePVT and TheGVT are declared in gvt.C, used everywhere
extern CkGroupID ThePVT;  
extern CkGroupID TheGVT;

class SRtable;  // from srtable.h
class SRentry;  // from srtable.h

/// Message to send info to GVT 
/** PVT sends processor virtual time and send/recv information to GVT.  
    GVT also sends info to next GVT index for next GVT invocation. */
class UpdateMsg : public CMessage_UpdateMsg {
public:
  /// PVT of local optimistic objects
  /** Used to send estimated GVT from one GVT invocation to next */
  POSE_TimeType optPVT;
  /// PVT of local conservative objects
  POSE_TimeType conPVT;
  /// Max timestamp in SR table
  POSE_TimeType maxSR;
  /// # sends/recvs at particular timestamps <= PVT
  SRentry *SRs;
  /// Count of entries in SRs
  int numEntries;
  /// Inactive status (GVT only)
  int inactive;
  /// Inactive time (GVT only)
  POSE_TimeType inactiveTime;
  /// Iterations of GVT since last LB (GVT only)
  int nextLB;
  /// Flag used by runGVT to call computeGVT indicating readiness
  int runGVTflag;
};

/// Message to send GVT estimate back to PVT
class GVTMsg : public CMessage_GVTMsg {
public:
  /// GVT estimate
  POSE_TimeType estGVT;
  /// Termination flag
  int done;
};

/// Prioritized int msg; used to force GVT calculations
class prioBcMsg : public CMessage_prioBcMsg {
public:
  int bc;
};

/// PVT chare group for computing local processor virtual time 
/** Keeps track of local sends/recvs and computes processor virtual time.
    Interacts with GVT to obtain new estimate and invokes fossil 
    collection and forward execution on objects with new estimate. */
class PVT : public Group {  
 private:
#ifdef POSE_STATS_ON
  localStat *localStats;
#endif
  /// PVT of local optimistic posers
  POSE_TimeType optPVT;
  /// PVT of local conservative posers
  POSE_TimeType conPVT;
  /// Last GVT estimate
  POSE_TimeType estGVT;       
  /// Last reported PVT estimate
  POSE_TimeType repPVT;       
  /// Simulation termination flag
  int simdone;
  /// Minimum send/recv timestamp in this iteration
  POSE_TimeType iterMin;
  /// Flag to indicate waiting for first send/recv of next iteration
  /** Used to indicate when to restructure the SendsAndRecvs table */
  int waitForFirst;
  /// Table to store send/recv timestamps
  SRtable *SendsAndRecvs;            
  /// List of objects registered with this PVT branch
  pvtObjects objs;           
  /// reduction-related vars
  int reportTo, reportsExpected, reportReduceTo, reportEnd;
  /// where the centralized GVT goes
  int gvtTurn;
  int specEventCount, eventCount;
  /// startPhase active flag
  int startPhaseActive;
 public:
  /// Basic Constructor
  PVT(void);
  PVT(CkMigrateMessage *) { };
  /// ENTRY: runs the PVT calculation and reports to GVT
  void startPhase(prioBcMsg *m);             
  /// ENTRY: receive GVT estimate; wake up objects
  /** Receives the new GVT estimate and termination flag; wakes up objects
      for fossil collection and forward execution with new GVT estimate. */
  void setGVT(GVTMsg *m);            
  /// Returns GVT estimate
  POSE_TimeType getGVT() { return estGVT; }    

  int getSpecEventCount() { return specEventCount; }    
  int getEventCount() { return eventCount; }    
  void incSpecEventCount() { specEventCount++; }    
  void incEventCount() { eventCount++; }
  void decEventCount() { eventCount--; }
  /// Returns termination flag
  int done() { return simdone; }
  /// Register poser with PVT
  int objRegister(int arrIdx, POSE_TimeType safeTime, int sync, sim *myPtr);
  /// Unregister poser from PVT
  void objRemove(int pvtIdx);
  /// Update send/recv table at timestamp
  void objUpdate(POSE_TimeType timestamp, int sr); 
  /// Update PVT with safeTime and send/recv table at timestamp
  void objUpdateOVT(int pvtIdx, POSE_TimeType safeTime, POSE_TimeType ovt);
  /// Reduction point for PVT reports
  void reportReduce(UpdateMsg *);
  /// Adds incoming send/recv information to a list
  void addSR(SRentry **SRs, SRentry *e, POSE_TimeType og, int ne);
  int getNumObjs() { return objs.getNumObjs(); }
};

/// GVT chare group for estimating GVT
/** Responsibility for GVT estimation shifts between branches after each
    GVT invocation. */
class GVT : public Group { 
private:
#ifdef POSE_STATS_ON
  localStat *localStats;
#endif
  /// Latest GVT estimate
  POSE_TimeType estGVT; 
  /// Inactivity status: number of iterations since GVT has changed
  int inactive;
  /// Time at which GVT last went inactive
  POSE_TimeType inactiveTime;
  /// Number of GVT iterations since last LB run
  int nextLBstart; 
  /// Earliest send/recv timestamp in previous GVT invocation
  POSE_TimeType lastEarliest;
  /// Number of sends at lastEarliest
  int lastSends;
  /// Number of receives at lastEarliest
  int lastRecvs;
  /// Number of PVT reports expected (1 or 2)
  int reportsExpected;
public:
  /// Basic Constructor
  GVT(void);
  GVT(CkMigrateMessage *) { };
  //Use this for Ccd calls
  //static void _runGVT(UpdateMsg *);
  /// ENTRY: Run the GVT
  /** Updates data fields with info from previous GVT estimation.  Fires
      off PVT calculations on all PVT branches. */
  void runGVT(UpdateMsg *);
  /// ENTRY: Gathers PVT reports; calculates and broadcasts GVT to PVTs
  void computeGVT(UpdateMsg *); 
  /// Adds incoming send/recv information to a list
  void addSR(SRentry **SRs, SRentry *e, POSE_TimeType og, int ne);
};

#endif
