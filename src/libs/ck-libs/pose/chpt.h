/// Checkpointing representation to be used with optimistic strategy
#ifndef CHPT_H
#define CHPT_H

/// Templated checkpointing class derived from rep
/** This class makes it possible for optimistic synchronization strategies
    to automatically checkpoint objects of classes derived from this type. */
template<class StateType> class chpt : public rep {
 public:
  int sinceLast;
  /// Basic Constructor
  chpt() { sinceLast = STORE_RATE; }
  /// Destructor
  virtual ~chpt() { }
  void registerTimestamp(int idx, eventMsg *m, int offset);
  /// Checkpoint the state
  void checkpoint(StateType *data);          
  /// Restore the state from a checkpoint
  /** Used during a rollback by the Undo method:  if the event being undone is
      not the final destination, we simply remove the checkpointed data; if the
      event is the final target for the rollback, then we restore the data, and
      remove the checkpointed data for it (it will be regenerated when the 
      target event gets re-executed). */
  void restore(StateType *data);
  virtual void pup(PUP::er &p) { rep::pup(p); }
  virtual void dump() { rep::dump(); }
};


/// Timestamps event message, sets priority, and records in spawned list
template<class StateType> 
void chpt<StateType>::registerTimestamp(int idx, eventMsg *m, int offset)
{
  m->Timestamp(ovt+offset);
  m->setPriority(ovt+offset-INT_MAX);
  parent->registerSent(ovt+offset);
  ((opt *)myStrat)->AddSpawnedEvent(idx, m->evID, m->timestamp);
}

/// Checkpoint the state
template<class StateType>
void chpt<StateType>::checkpoint(StateType *data)
{
#ifdef POSE_STATS_ON
  localStat *localStats = (localStat *)CkLocalBranch(theLocalStats);
  localStats->SwitchTimer(CP_TIMER);
#endif
  CmiAssert(!(parent->myStrat->currentEvent->cpData));
  //  if ((myStrat->currentEvent->timestamp > 
  //myStrat->currentEvent->prev->timestamp) || (sinceLast == STORE_RATE)) {
  if ((sinceLast == STORE_RATE) || (CpvAccess(stateRecovery) == 1) ||
      (myStrat->currentEvent->prev == parent->eq->front())) {
    myStrat->currentEvent->cpData = new StateType;
    myStrat->currentEvent->cpData->copy = 1;
    *((StateType *)myStrat->currentEvent->cpData) = *data;
    sinceLast = 0;
  }
  else sinceLast++;
#ifdef POSE_STATS_ON
  localStats->SwitchTimer(SIM_TIMER);
#endif
}

/// Restore the state from a checkpoint
template<class StateType> 
void chpt<StateType>::restore(StateType *data) 
{
  if (myStrat->currentEvent == myStrat->targetEvent) {
    if (myStrat->targetEvent->cpData) {
      *data = *((StateType *)myStrat->targetEvent->cpData);
      delete myStrat->targetEvent->cpData;
      myStrat->targetEvent->cpData = NULL;
    }
  }
  if (myStrat->currentEvent->cpData) {
    delete myStrat->currentEvent->cpData;
    myStrat->currentEvent->cpData = NULL;
  }
}

#endif
