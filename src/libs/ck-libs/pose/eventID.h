/// Unique event IDs for POSE events
/** Provides the event ID data structure and function GetEventID to generate
    the unique event IDs */
#ifndef EVENTID_H
#define EVENTID_H
#include "charm++.h"
#include "limits.h"

#ifdef WIN32
#define snprintf _snprintf
#endif

/// Unique identifier for a POSE event
class eventID 
{
  /// PE identifier field ensures uniqueness across PEs
  int pe;
  /// Control field for ordering events with same tiemstamp
  int control;
 public:
  /// Large number for identifier unique on creation PE
  unsigned int id; 
  /// Basic Constructor
  eventID() { id = 0; pe = CkMyPe(); control = 0; }          
  /// Init to default unused values
  void init() { id = 0; pe = -1; control = INT_MIN; }
  /// Get next value for eventID
  /** increments id field for this eventID */
  void incEventID() {
    id++;
    if (id == 0) CkPrintf("WARNING: event ID rollover occurred!\n");
  }
  /// Assignment operator
  eventID& operator=(const eventID& e) { 
    CmiAssert((e.pe >= 0) || (e.pe < CkNumPes()));
    id = e.id;  pe = e.pe;  control = e.control; return *this;
  }
  /// get source PE
  int getPE() { return pe; }
  /// get source control
  int getControl() { return control; }
  /// set source control
  void setControl(int ctrl) { control = ctrl; }
  /// Equality comparison operator
  int operator==(const eventID& o) {
    return ((id == o.id) && (pe == o.pe) && (control == o.control));
  }
  /// Less than/equality comparison operator
  /** Provides a way to sort events by event ID */
  int operator<=(const eventID& o) {
    return (control <= o.control);
  }
  /// Less than comparison operator
  /** Provides a way to sort events by event ID */
  int operator<(const eventID& o) {
    return (control < o.control);
  }
  /// Greater than/equality comparison operator
  /** Provides a way to sort events by event ID */
  int operator>=(const eventID& o) {
    return (control >= o.control);
  }
  /// Greater than comparison operator
  /** Provides a way to sort events by event ID */
  int operator> (const eventID& o) {
    return (control >  o.control);
  }

  /// Dump all data fields
  void dump() { 
    CmiAssert((pe >= 0) && (pe < CkNumPes())); 
    CkPrintf("%d.%d", id, pe); 
  }    
  char *sdump(char *s) { sprintf(s, "%d.%d", id, pe); return s;}    
  char *sndump(char *s,size_t n) { snprintf(s,n,"%d.%d", id, pe); return s;}
  /// Pack/unpack/sizing operator
  void pup(class PUP::er &p) { p(id); p(pe); }  
  void sanitize() { CkAssert((pe > -1) && (pe < CkNumPes())); }
};

/// Generates and returns unique event IDs
const eventID& GetEventID();                    

#endif
