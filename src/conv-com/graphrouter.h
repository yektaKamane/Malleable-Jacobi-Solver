/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

#ifndef _GRAPHROUTER_H
#define _GRAPHROUTER_H

#include <math.h>
#include <converse.h>
#include "comlib.h"
#include "petable.h"

#include "persistent.h"

class TopologyDescriptor {
 public:
    //TopologyDescriptor(int npes, int mype)=0;
    //Entry methods which will define the specific graph.
    virtual void getNeighbors(int &npes, int *pelist) = 0;
    virtual int getNumStages() = 0;
    virtual int getNumSteps(int stage) = 0;
    virtual void getPesToSend(int step, int stage, int &npesToSend, 
                              int *pelist, int &nextpe) = 0;
    virtual int getNumMessagesExpected(int stage) = 0;
};

//A generalized virtual topology based router. To make it a specific
//router a new topology class needs to br created and passed to this.
class GraphRouter : public Router
{
    PeTable *PeGraph;
    TopologyDescriptor *tp;
    
    int *pesToSend, *gpes;
    int numNeighbors, *neighborPeList, nstages;
    
    int MyPe, NumPes, currentIteration;
    int *recvExpected, *recvCount;
    int curStage, *stageComplete;

    void sendMessages(comID id, int stage);
    void init(int numPes, int myPe, TopologyDescriptor *tp);
    
#if CMK_PERSISTENT_COMM
    PersistentHandle *handlerArrayOdd, *handlerArrayEven;
#endif          
    
 public:
    GraphRouter(int numPes, int myPe);
    //GraphRouter(int numPes, int myPe, int topid);
    ~GraphRouter();
  
    //Router enrty methods overridden here.

    void NumDeposits(comID, int);
    void EachToAllMulticast(comID , int , void *, int);
    void EachToManyMulticast(comID , int , void *, int, int *, int);
    void RecvManyMsg(comID, char *);
    void ProcManyMsg(comID, char *);
    void DummyEP(comID id, int);
    
    void SetMap(int *);
};

#endif
