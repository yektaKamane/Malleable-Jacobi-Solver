//Contains the structures for bg logs and timelines
#ifndef BLUE_LOGS_H
#define BLUE_LOGS_H

#include <string.h>

#include "blue_defs.h"
#include "cklists.h"

extern int bgcorroff;

// identifier for a message which records the source node that generate
// this message and a message sequence number (_msgID)
class BgMsgID
{
private:
  int _node;		// node number where the message is created
  int _msgID;		// local index number on pe

public:
  BgMsgID(): _node(-1), _msgID(-1) {}
  BgMsgID(int p, int m): _node(p), _msgID(m) {}
  void pup(PUP::er &p) {
    p|_node; p|_msgID;
  }
  inline int operator == (const BgMsgID &m) {
    return _node == m._node && _msgID == m._msgID;
  }
  int node() { return _node; }
  int msgID() { return _msgID; }
};

/**
  a message sent event in timeline
*/
class BgMsgEntry {
  friend class BgTimeLog;
public:
  int msgID;
  int dstPe;          // dest bg node in global sequence
  double sendTime;	// msg sending offset in the event
  double recvTime;	// predicted recv time with delay
#if DELAY_SEND
  char *sendMsg;	// real msg
#endif
  CmiInt2 tID;		// destination worker thread ID
  int msgsize;		// message size
private:
  BgMsgEntry() {}
public:
  BgMsgEntry(char *msg, int node, int tid, int local);
  inline void print() {
    CmiPrintf("msgID:%d sent:%f recvtime:%f dstPe:%d\n", msgID, sendTime, recvTime, dstPe);
  }
  void write(FILE *fp) {
    fprintf(fp, "msgID:%d sent:%f recvtime:%f dstPe:%d\n", msgID, sendTime, recvTime, dstPe);
  }
#if DELAY_SEND
  void send();
#endif
  void pup(PUP::er &p) {
    p|msgID; p|dstPe; p|sendTime; p|recvTime; p|tID; p|msgsize;
  }
};


/**
  event for higher level of tracing like trace projections
*/
class bgEvents {
private:
  bgEventCallBackFn  callbackFn;
  void* usrPtr;
public:
  void*   data;         // e.g. can be pointer to trace projection log entry
  int     index;		// index of the event to its original log pool.
  double  rTime;	// relative time from the start entry
  char   eType;
  bgEvents(): index(-1) {}
  bgEvents(void *d, int idx, double t, bgEventCallBackFn fn, void *ptr, char e):
	data(d), index(idx), rTime(t), callbackFn(fn), usrPtr(ptr), eType(e) {}
  inline void update(double startT, double recvT, int e) {
	if (eType==e) callbackFn(data, startT+rTime, recvT, usrPtr);
  }
  void print();
  void write(FILE *fp);
  void pup(PUP::er &p);
};

#define BG_STARTSIM     0x1

class BgTimeLineRec;
/**
  one time log for an handler function;
  it record a list of message sent events in an execution of handler
*/
class BgTimeLog {
public:
  int ep;
  int seqno;
  BgMsgID  msgId;	// incoming message that generates this log

  double recvTime;	//Time at which the message was received in 'inbuffer'
  double startTime, endTime;
  double oldStartTime, execTime;
  double effRecvTime;

//  int threadNum;	// by guna, for seq load balancing  ???

  CkVec< BgMsgEntry * > msgs;
  CkVec< bgEvents * > evts;
  CkVec< BgTimeLog* > backwardDeps;
  CkVec< BgTimeLog* > forwardDeps;
  char doCorrect;
  char flag;
  char name[20];

  friend class BgTimeLineRec;
public:
  BgTimeLog(BgTimeLog *);
  BgTimeLog(char *msg, char *str=NULL);
  BgTimeLog(): ep(-1), recvTime(.0), startTime(.0), endTime(.0),
	       execTime(.0), effRecvTime(INVALIDTIME), seqno(0), doCorrect(1) 
    {strcpy(name,"dummyname");}
  BgTimeLog(int epc, char* name, double sTime, double eTime);
  BgTimeLog(int epc, char* name, double sTime);
  ~BgTimeLog();

  inline void setExecTime() {
           execTime = endTime - startTime;
           if(execTime < BG_EPSILON && execTime > -BG_EPSILON)
             execTime = 0.0;
           CmiAssert(execTime >= 0.0);
         }
  inline void addMsg(char *msg, int node, int tid, int local) { 
           msgs.push_back(new BgMsgEntry(msg, node, tid, local)); 
         }
  void closeLog();
  void print(int node, int th);
  void write(FILE *fp);

  inline void setStartEvent() { flag |= BG_STARTSIM; }
  inline int isStartEvent() { return (flag & BG_STARTSIM); }

  // add backward dep of the log corresponent to msg
  void addMsgBackwardDep(BgTimeLineRec &tlinerec, void* msg);
  void addBackwardDep(BgTimeLog* log);
  //takes a list of Logs on which this log is dependent (backwardDeps) 
  void addBackwardDeps(CkVec<BgTimeLog*> logs);
  void addBackwardDeps(CkVec<void*> logs);
  int bDepExists(BgTimeLog* log);			// by guna
  //Returns earliest time by which all backward dependents ended  
  // return the last eff recv time
  double getEndOfBackwardDeps() {
    double maxEndTime =0.0;
    for(int i=0;i<backwardDeps.length();i++)
//    maxEndTime = max(maxEndTime,backwardDeps[i]->endTime);
      maxEndTime = BG_MAX(maxEndTime,backwardDeps[i]->effRecvTime);
      
    return maxEndTime;
  }

  inline void addEvent(void *data,int idx,double absT,bgEventCallBackFn fn,void *p,int e) { 
    evts.push_back(new bgEvents(data, idx, absT-startTime, fn, p, e)); 
  }
  inline void updateEvents(int e) {
    for (int i=0; i<evts.length(); i++)
      evts[i]->update(startTime ,recvTime, e);
  }
  inline double key() { return effRecvTime; }
  inline int compareKey(BgTimeLog* otherLog){
    if(((isZero(effRecvTime-otherLog->effRecvTime))&&(seqno < otherLog->seqno))
       ||(isLess(effRecvTime,otherLog->effRecvTime)))
      return -1;
    return 1;
  }
  inline int isEqual(BgTimeLog* otherLog){
    return (otherLog==this);
  }
  void pup(PUP::er &p);

#if DELAY_SEND
  void send() {
    for (int i=0; i<msgs.length(); i++)
      msgs[i]->send();
  }
#endif
};


/**
  Timeline for a VP
*/
typedef CkQ< BgTimeLog *> BgTimeLine;

/**
  A wrapper for CkQ of BgTimeLine
*/
class BgTimeLineRec {
public:
  BgTimeLine  timeline;
  int         commit;
  int         startIdx;
  int         startCorrFlag;
  int         correctSendIdx;
  int 	      counter;
  double      minCorrection;
  BgTimeLog  *bgCurLog;		/* current unfinished log */
  BgTimeLog  *bgPrevLog;	/* previous log that should make dependency */
#if DELAY_SEND
  CkQ<BgTimeLog *>   sendingLogs;	// send buffered
#endif
public:
  BgTimeLineRec(): timeline(1024), commit(0), counter(1), correctSendIdx(0), 
		   startIdx(0), bgCurLog(NULL), bgPrevLog(NULL) {
      if (bgcorroff) startCorrFlag=0; else startCorrFlag=1;
      minCorrection = INVALIDTIME;
    }
  ~BgTimeLineRec() {
      for (int i=0; i<timeline.length(); i++)  delete timeline[i];
    }
  BgTimeLog * operator[](size_t n) {
	CmiAssert(n!=(size_t)-1);
        return timeline[n];
    }
  int length() { return timeline.length(); }
  // special enq which will assign seqno
  void enq(BgTimeLog *log, int isnew) {
	log->seqno = counter++;
  	timeline.enq(log);
#if DELAY_SEND
	if (isnew) sendingLogs.enq(log);
#endif
  }
  void setStartIdx(){
    startIdx = timeline.length();
  }
  double computeUtil(int *numRealMsgs){
    //From startIdx to the end of the timeline
    double total=0.0;
    int tlineLen = length();
    for(int i=0;i<tlineLen;i++) {
      BgTimeLog *log = timeline[i];
      total += log->execTime;
      *numRealMsgs += log->msgs.length();
    }
    return total;
  }
  inline void clearSendingLogs() {
#if DELAY_SEND
    while (!sendingLogs.isEmpty()) {
      BgTimeLog *log = sendingLogs.deq();
      log->send();
    }
#endif
  }
  void logEntryStart(char *m);
//  void logEntryCommit();
  void logEntryInsert(BgTimeLog* log);
  void logEntryStart(BgTimeLog* log);
  void logEntryClose();
  void logEntrySplit();
  BgTimeLog *getTimeLogOnThread(const BgMsgID &msgId, int *index);

  void pup(PUP::er &p);
};

// BigSim log function API
int BgIsInALog(BgTimeLineRec &tlinerec);
BgTimeLog *BgLastLog(BgTimeLineRec &tlinerec);
void BgAddBackwardDep(BgTimeLog *curlog, BgTimeLog* deplog);
BgTimeLog *BgStartLogByName(BgTimeLineRec &tlinerec, int ep, char *name, double starttime, BgTimeLog *prevLog);
void BgEndLastLog(BgTimeLineRec &tlinerec);

int BgLoadTraceSummary(char *fname, int &totalProcs, int &numX, int &numY, int &numZ, int &numCth, int &numWth, int &numPes);
void BgReadProc(int procNum, int numWth ,int numPes, int totalProcs, int* allNodeOffsets, BgTimeLineRec& tlinerec);
int* BgLoadOffsets(int totalProcs, int numPes);
void BgWriteThreadTimeLine(char *fname, int x, int y, int z, int th, BgTimeLine &tline);

#endif
