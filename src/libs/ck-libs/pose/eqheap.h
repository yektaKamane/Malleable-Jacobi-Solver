/// Heap structure used for unexecuted portion of the eventQueue
/** This should provide rapid insertion/deletion of new events to the 
    event queue */
#ifndef EQHEAP_H
#define EQHEAP_H

/// Structure for storing events on a heap
class HeapNode
{
 public:
  /// Size of subheap topped by this node
  int subheapsize;         
  /// The event stored here
  Event *e;                
  /// Left and right subheaps
  HeapNode *left, *right;  
  /// Basic Constructor
  HeapNode() { subheapsize = 0; e = NULL; left = right = NULL; }
  /// Initializing Constructor
  HeapNode(Event *ev, int sz, HeapNode *l, HeapNode *r) { 
    subheapsize = sz; e = ev; left = l; right = r; 
  }
  /// Destructor
  ~HeapNode(){if (left) delete left; if (right) delete right; if (e) delete e;}
  /// Insert event in heap
  /** Insert event e in this subheap; designed to find insertion position
      quickly, at the cost of creating unbalanced or high & narrow heaps */
  void insert(Event *e);                    
  /// Join this heap with h
  /** Join this heap with h and return the new heap; uses quickest join method 
      possible at expense of creating unbalanced tree */
  HeapNode *conjoin(HeapNode *h);           
  /// Remove heap node matching evID
  int remove(eventID evID, POSE_TimeType timestamp);  
  /// Dump all data fields
  void dump();                 
  /// Pack/unpack/sizing operator
  /** Packs/sizes the entire heap, DOES NOT UNPACK HEAP!!! */
  void pup(PUP::er &p);                     
};

/// Heap structure to store events in unexecuted portion of event queue
class EqHeap {  
  /// Size of heap
  int heapSize; 
 public:
  /// Top node of heap  
  HeapNode *top;
  /// Basic Constructor
  EqHeap() { heapSize = 0; top = NULL; }
  /// Destructor
  ~EqHeap() { if (top) delete top; } 
  /// Insert event e in heap with low timestamps at top of heap
  void InsertEvent(Event *e);              
  /// Return event on top of heap, deleting it from the heap
  /** Returns event at top of heap if one exists, null otherwise; deletes top
      node in heap, conjoining left and right subheaps */
  Event *GetAndRemoveTopEvent();           
  /// Delete event from heap
  /** Delete the node with event corresponding to evID and timestamp; 
      returns 1 if an event was successfully deleted, 0 if the event was not
      found in the heap */
  int DeleteEvent(eventID evID, POSE_TimeType timestamp);  
  /// Dump entire heap
  void dump();                      
  /// Pack/unpack/sizing operator
  /** Pups entire heap relying on recursive HeapNode::pup */
  void pup(PUP::er &p);     
  /// Check validity of data fields
  void sanitize();
};

#endif
