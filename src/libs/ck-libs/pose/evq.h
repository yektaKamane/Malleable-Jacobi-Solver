/// Queue of executed and unexecuted events on a poser
#ifndef EVQ_H
#define EVQ_H

//#define EQ_SANITIZE 1

/// The event queue
/** Doubly-linked list with front and back sentinels and a heap of unexecuted
    events */
class eventQueue {
  /// This helper method cleans up all the commit code by removing the stats
  void CommitStatsHelper(Event *commitPtr);
 public:
  /// Sentinel nodes
  Event *frontPtr, *backPtr;
  /// Pointer to next unexecuted event
  Event *currentPtr;
  /// Event to rollback to
  Event *RBevent;
  /// Heap of unexecuted events
  EqHeap *eqh;
  /// Largest unexecuted event timestamp in queue
  POSE_TimeType largest;
  /// number of unexecuted events in the queue
  unsigned int eventCount;
  /// Output file name for stats for DOP calculation
  char filename[20];
  /// Output file pointer for stats for DOP calculation
  FILE *fp;
  /// Coarse memory usage
  unsigned int mem_usage;
  /// Keep track of last logged VT for this object so no duplicates are logged
  POSE_TimeType lastLoggedVT;
  /// Basic Constructor
  /** Creates front and back sentinel nodes, connects them, and inits pointers
      and heap. */
  eventQueue();
  /// Destructor
  ~eventQueue();
  /// Insert e in the queue in timestamp order
  /** If executed events with same timestamp exist, insert e at the back of 
      these, returns 0 if rollback necessary, 1 otherwise. */
  void InsertEvent(Event *e);      
  /// Insert e in the queue in timestamp order
  /** If executed events with same timestamp exist, sort e into them based on
      eventID, returns 0 if rollback necessary, 1 otherwise. */
  void InsertEventDeterministic(Event *e);      
  /// Return front pointer
  Event *front() { return frontPtr; }
  /// Return back pointer
  Event *back() { return backPtr; }
  /// Move currentPtr to next event in queue
  /** If no more events, take one from heap */
  void ShiftEvent() { 
    Event *e;
#ifdef EQ_SANITIZE
    sanitize();
#endif
    CmiAssert(currentPtr->next != NULL);
    currentPtr = currentPtr->next; // set currentPtr to next event
    if ((currentPtr == backPtr) && (eqh->top)) { // currentPtr on back sentinel
      e = eqh->GetAndRemoveTopEvent(); // get next event from heap
      // insert event in list
      e->prev = currentPtr->prev;
      e->next = currentPtr;
      currentPtr->prev = e;
      e->prev->next = e;
      currentPtr = e;
    }
    if (currentPtr == backPtr) largest = POSE_UnsetTS;
    else FindLargest();
    eventCount--;
#ifdef EQ_SANITIZE
    sanitize();
#endif
  }               
  /// Commit (delete) events before target timestamp ts
  void CommitEvents(sim *obj, POSE_TimeType ts); 
  /// Commit (delete) all events
  void CommitAll(sim *obj); 
  /// Change currentPtr to point to event e
  /** Be very very careful with this -- avoid using if possible */
  void SetCurrentPtr(Event *e);
  /// Delete event and reconnect surrounding events in queue
  void DeleteEvent(Event *ev);
  /// Return the event ID of the event pointed to by currentPtr
  const eventID& CurrentEventID() { return currentPtr->evID; }
  /// Add id, e and ts as an entry in currentPtr's spawned list
  /** The poser e was sent to is id */
  void AddSpawnToCurrent(int id, eventID e, POSE_TimeType ts) {
    SpawnedEvent *newnode = 
      new SpawnedEvent(id, e, ts, currentPtr->spawnedList);
    CmiAssert(currentPtr->done == 2);
    currentPtr->spawnedList = newnode;
  }
  /// Return the first entry in currentPtr's spawned list and remove it
  SpawnedEvent *GetNextCurrentSpawn() {
    SpawnedEvent *tmp = currentPtr->spawnedList;
    if (tmp) currentPtr->spawnedList = tmp->next;
    return tmp;
  }
  /// Find the largest timestamp of the unexecuted events
  void FindLargest() {
    POSE_TimeType hs = eqh->FindMax();
    if (backPtr->prev->done == 0) largest = backPtr->prev->timestamp;
    else largest = POSE_UnsetTS;
    if (largest < hs) largest = hs;
  }
  /// Set rollback point to event e
  void SetRBevent(Event *e) {
    if (!RBevent) RBevent = e; 
#ifdef DETERMINISTIC_EVENTS
    else if ((RBevent->timestamp > e->timestamp) ||
	     ((RBevent->timestamp == e->timestamp) && 
	      (e->evID < RBevent->evID))) {
      CmiAssert(RBevent->prev->next == RBevent);
      CmiAssert(RBevent->next->prev == RBevent);
      RBevent = e;
    }
#else
    else if (RBevent->timestamp > e->timestamp) {
      CmiAssert(RBevent->prev->next == RBevent);
      CmiAssert(RBevent->next->prev == RBevent);
      RBevent = e;
    }
#endif
  }
  /// Dump the event queue
  void dump();        
  /// Pack/unpack/sizing operator
  void pup(PUP::er &p);   
  /// Check validity of data fields
  void sanitize();
};

#endif
