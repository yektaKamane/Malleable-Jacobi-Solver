/// SendRecvTable for POSE GVT calculations
#ifndef SRTABLE_H
#define SRTABLE_H
#include "pose.h"

#define MAX_B 10

//#define SR_SANITIZE 1

class UpdateMsg; // from gvt.h

/// An entry for storing the number of sends and recvs at a timestamp
/** This class is used in SendRecvTable to store residual send/recv data */
class SRentry {
 public:
  /// Timestamp of the message
  POSE_TimeType timestamp;  
  /// The number of messages sent with theTimestamp
  int sends;
  /// The number of messages sent with theTimestamp
  int recvs;
  /// Next SRentry in list
  SRentry *next;
  /// Basic constructor
  /** Initializes all data members */
  SRentry() { 
    timestamp = POSE_UnsetTS; sends = recvs = 0; next = NULL; 
  }
  /// Initializing constructor 1
  /** Initializes timestamp & next w/parameters, sends & recvs to 0 */
  SRentry(POSE_TimeType ts, SRentry *p) { 
    timestamp = ts; sends = recvs = 0; next = p; 
  }
  /// Initializing constructor 2
  /** Initializes timestamp, send/recv count and next
      w/parameters */
  SRentry(POSE_TimeType ts, int sr, SRentry *p) {
    timestamp = ts;  next = p; 
    if (sr == SEND) { sends = 1; recvs = 0; }
    else { sends = 0; recvs = 1; }
  }
  /// Initializing constructor
  /** Initializes timestamp and send/recv count w/parameters */
  SRentry(POSE_TimeType ts, int sr) {
    timestamp = ts;  next = NULL; 
    if (sr == SEND) { sends = 1; recvs = 0; }
    else { sends = 0; recvs = 1; }
  }
  /// Assignment operator
  SRentry& operator=(const SRentry& e) {
    timestamp = e.timestamp;
    sends = e.sends;
    recvs = e.recvs;
    return *this;
  }
  /// Dump data fields
  void dump() {
    if (next)
      CkPrintf("TS:%d #s:%d #r:%d n:!NULL ", timestamp, sends, recvs); 
    else CkPrintf("TS:%d #s:%d #r:%d n:NULL ",timestamp, sends, recvs);
  }
  /// Check validity of data fields
  void sanitize() {
    CmiAssert(timestamp >= POSE_UnsetTS); // should be POSE_UnsetTS or > if initialized
    CmiAssert(sends >= 0);  // cannot be less than zero
    CmiAssert(recvs >= 0);  // cannot be less than zero
    if (next == NULL) return;   // next can be NULL
    // if next != NULL, check if pointer looks valid
    POSE_TimeType test_ts = next->timestamp;
    int test_sendCount = next->sends;
    int test_recvCount = next->recvs;
    SRentry *test_next = next->next;
    next->timestamp = test_ts;
    next->sends = test_sendCount;
    next->recvs = test_recvCount;
    next->next = test_next;
  }
};

/// A table for storing the number of sends and recvs at a timestamp
/** This class is used in GVT to keep track of messages sent/received */
class SRtable {
 public:
  /// Base timestamp to index tables
  /** offset is the current GVT */
  POSE_TimeType offset;
  /// Number of buckets to sort sends/recvs into
  /** Recomputed with each new offset */
  int b;
  /// Size of each bucket
  /** Recomputed with each new offset */
  int size_b;
  /// The buckets to sort sends/recvs into
  /** Only entries [0..b-1] are used */
  SRentry *buckets[MAX_B];
  /// Pointers to the last entry of each bucket
  /** Only entries [0..b-1] are used */
  SRentry *end_bucket[MAX_B];
  /// Error checking on bucket counts
  int sends[MAX_B], recvs[MAX_B], ofSends, ofRecvs;
  /// The overflow bucket
  /** What doesn't fit in buckets goes here */
  SRentry *overflow;
  /// End entry of overflow
  SRentry *end_overflow;
  /// Number of distinct timestamp entries per bucket
  /** This is computed in CompressAndSortBucket */
  int numEntries[MAX_B];
  /// Number of distinct entries in overflow bucket
  /** This is computed in CompressAndSortBucket */
  int numOverflow;
  
  /// Basic constructor
  /** Initializes all data fields */
  SRtable();
  /// Destructor
  ~SRtable() { FreeTable(); }
  /// Initialize table to a minimum size
  void Initialize();
  /// Insert send/recv record sr at timestamp ts
  void Insert(POSE_TimeType ts, int sr) {
#ifdef SR_SANITIZE
    sanitize();
#endif
    CmiAssert(ts >= offset);
    CmiAssert((sr == 0) || (sr == 1));
    if (size_b == -1) size_b = 1 + ts/b;
    POSE_TimeType destBkt = (ts-offset)/size_b;  // which bucket?
    SRentry *e = new SRentry(ts, sr, NULL);
    if (destBkt >= b) { // put in overflow bucket
      if (overflow) {
	if (end_overflow->timestamp == ts) { // an entry at ts exists
	  if (sr == SEND) end_overflow->sends++;
	  else end_overflow->recvs++;
	  delete e;
	}
	else { // no entry with that timestamp is handy
	  end_overflow->next = e;
	  end_overflow = e;
	}
      }
      else overflow = end_overflow = e;
      if (sr == SEND) ofSends++;
      else ofRecvs++;
    }
    else { // put in buckets[destBkt]
      if (buckets[destBkt]) {
	if (end_bucket[destBkt]->timestamp == ts) { 
	  // an entry at that timestamp exists
	  if (sr == SEND) end_bucket[destBkt]->sends++;
	  else end_bucket[destBkt]->recvs++;
	  delete e;
	}
	else { // no entry with that timestamp is handy
	  end_bucket[destBkt]->next = e;
	  end_bucket[destBkt] = e;
	}
      }
      else buckets[destBkt] = end_bucket[destBkt] = e;
      if (sr == SEND) sends[destBkt]++;
      else recvs[destBkt]++;
    }
#ifdef SR_SANITIZE
    sanitize();
#endif
  } 
  /// Insert an existing SRentry e
  void Insert(SRentry *e) {
#ifdef SR_SANITIZE
    sanitize();
#endif
    CmiAssert(e != NULL);
    CmiAssert(e->timestamp >= offset);
    POSE_TimeType destBkt = (e->timestamp-offset)/size_b;
    e->next = NULL;
    if (destBkt >= b) { // put in overflow bucket
      if (overflow) {
	end_overflow->next = e;
	end_overflow = e;
      }
      else overflow = end_overflow = e;
      ofSends += e->sends;
      ofRecvs += e->recvs;
    }
    else { // put in buckets[destBkt]
      if (buckets[destBkt]) {
	end_bucket[destBkt]->next = e;
	end_bucket[destBkt] = e;
      }
      else buckets[destBkt] = end_bucket[destBkt] = e;
      sends[destBkt] += e->sends;
      recvs[destBkt] += e->recvs;
    }
#ifdef SR_SANITIZE
    sanitize();
#endif
  }
  /// Restructure the table according to new GVT estimate and first send/recv
  /** Number of buckets and bucket size are determined from firstTS, and
      entries below newGVTest are discarded. */
  void Restructure(POSE_TimeType newGVTest, POSE_TimeType firstTS,int firstSR);
  /// Move contents of bkt to new bucket structure
  void MapToBuckets(SRentry *bkt, SRentry *endBkt, int *s, int *r);
  /// Compress and pack table into an UpdateMsg and return it
  UpdateMsg *PackTable(POSE_TimeType pvt);
  /// CompressAndSort all buckets
  void SortTable();
  /// Compress a bucket so all SRentries have unique timestamps and are sorted
  void CompressAndSortBucket(POSE_TimeType i, int is_overflow);
  /// Free all buckets and overflows, reset all counts
  void FreeTable();
  /// Dump data fields
  void dump();    
  /// Check validity of data fields
  void sanitize();
  /// Test this class
  void self_test();
};

#endif
