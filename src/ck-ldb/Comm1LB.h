/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _Comm1LB_H_
#define _Comm1LB_H_

#include <CentralLB.h>
#include "Comm1LB.decl.h"

#include "CommLBHeap.h"
#include "CommLB.h"

#define CUT_OFF_FACTOR 1.200

struct alloc_struct{
  double load;
  int nbyte;
  int nmsg;
};

void CreateComm1LB();

class Comm1LB : public CentralLB {
public:
  int nobj,npe;
  alloc_struct ** alloc_array;
  graph * object_graph;
  Comm1LB();
  Comm1LB(CkMigrateMessage *m):CentralLB(m) {}
private:
  CmiBool QueryBalanceNow(int step);
  LBMigrateMsg* Strategy(CentralLB::LDStats* stats, int count);
  void alloc(int pe, int id, double load, int nmsg, int nbyte);
  double compute_cost(int id, int pe, int n_alloc, int &out_msg, int &out_byte); 
  void add_graph(int x, int y, int data, int nmsg);
};

#endif


/*@}*/

