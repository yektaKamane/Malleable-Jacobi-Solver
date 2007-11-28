/// Global POSE data and functions; includes and dependencies handled here
/** This code provides all the major global data structures plus a 
    coordination entity that handles initialization and termination of the
    simulation. */
#ifndef POSE_H
#define POSE_H


#include "pose.decl.h"
#include "pose_config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "charm++.h"
#include "eventID.h"
#include "mempool.h"
#include "memory_temporal.h"
#include "srtable.h"
#include "stats.h"
#include "cancel.h"

class eventMsg; // defined later in sim.h
class rep; // defined later in rep.h
#include "event.h"
#include "eqheap.h"

class sim; // defined later in sim.h
#include "pvtobj.h"
#include "lbObject.h"
#include "ldbal.h"
#include "gvt.h"
#include "evq.h"

class strat; // defined later in strat.h
#include "rep.h"
#include "strat.h"
#include "sim.h"
#include "sim.decl.h"
#include "opt.h"
#include "opt2.h"
#include "opt3.h"
#include "spec.h"
#include "adapt.h"
#include "adapt2.h"
#include "adapt3.h"
#include "adapt4.h"
#include "cons.h"
#include "seq.h"
#include "chpt.h"

/// Main initialization for all of POSE
void POSE_init(); 
void POSE_init(int ET); 
void POSE_init(int IDflag, int ET); 
/// Start POSE simulation timer
void POSE_startTimer(); 
/// Use Inactivity Detection to terminate program
void POSE_useID();

/// Use a user-specified end time to terminate program
/** Also uses inactivity detection in conjunction with end time */
void POSE_useET(POSE_TimeType et); 
/// Specify an optional callback to be called when simulation terminates
void POSE_registerCallBack(CkCallback cb);
/// Stop POSE simulation
/** Stops timer so statistics collection, callback, final output, etc. 
    are not counted in simulation time. */
void POSE_stop(); 
// Prepare to exit
void POSE_prepExit(void *param, void *msg); 
/// Exit simulation program
void POSE_exit(); 

/// User specified busy wait time (for grainsize testing)
extern double busyWait;

/// Simulation start time
extern double sim_timer;

/// Simulation end time
extern POSE_TimeType POSE_endtime;

/// Inactivity detection flag
extern int POSE_inactDetect;

/// Global clock (for sequential simulation)
extern POSE_TimeType POSE_GlobalClock;
extern POSE_TimeType POSE_GlobalTS;

/// For getting access to the commlib strategy
extern ComlibInstanceHandle POSE_commlib_insthndl;

extern POSE_Config pose_config;

/// Set busy wait time
void POSE_set_busy_wait(double n);
/// Busy wait for busyWait
void POSE_busy_wait();
/// Busy wait for n
void POSE_busy_wait(double n);

/// Flag to indicate how foward execution should proceed
/** 0 for normal forward execution; 1 for state recovery only */
CpvExtern(int, stateRecovery);

#ifdef POSE_COMM_ON
extern int comm_debug;
#endif

/// Class for user-specified callback
class callBack : public CMessage_callBack
{
 public:
  CkCallback callback;
};

/// Coordinator of simulation initialization, start and termination
class pose : public Chare {
 private:
  /// A callback to execute on termination
  /** If this is used, control is turned over to this at the very end of the
      simulation. */
  CkCallback cb;
  /// Flag to indicate if a callback will be used
  int callBackSet;
 public:
  /// Basic Constructor
  pose(void) : callBackSet(0) { 
#ifdef VERBOSE_DEBUG
  CkPrintf("[%d] constructing pose\n",CkMyPe());
#endif

 }
  pose(CkMigrateMessage *) { }
  /// Register the callback with POSE
  void registerCallBack(callBack *);
  /// Stop the simulation
  /** Stops timer and gathers POSE statistics and proceeds to exit. */
  void stop();
  /// Exit the simulation
  /** Executes callback before terminating program. */

  //! handle stats output before exiting if necessary 
  void prepExit();

  void exit();
};

#endif
