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

#ifndef _ORBLB_H_
#define _ORBLB_H_

#include "CentralLB.h"
#include "OrbLB.decl.h"

void CreateOrbLB();

class OrbLB : public CentralLB {
public:/* <- Sun CC demands Partition be public for ComputeLoad to access it. */

  class Partition {
  public:
    int refno;
    double load;			// total load in this set
    int origin[3];			// box coordinates
    int corner[3];
    int  count;
    int node, mapped;
  public:
    Partition(): refno(0), load(0.0), node(-1), mapped(0) {};
  };

private:  
  struct ComputeLoad {
    int id;
    int v[3];
    double load;
    int  refno;
    double  tv;
    Partition * partition;
  };
  
  
  struct VecArray {
    int v;
    int id;
  };
  
  enum {XDIR=0, YDIR, ZDIR};
  
  int P;
  ComputeLoad *computeLoad;
  int nObjs;
  VecArray  *(vArray[3]);
  Partition *partitions;
  Partition top_partition;
  int npartition;
  int currentp, refno;
  
  void strategy();
  void rec_divide(int, Partition&);
  void setVal(int x, int y, int z);
  int sort_partition(int x, int p, int r);
  void qsort(int x, int p, int r);
  void quicksort(int x);
  void mapPartitionsToNodes();

public:
  double overLoad;

public:
  OrbLB(const CkLBOptions &);
  OrbLB(CkMigrateMessage *m):CentralLB(m) {}
private:
  CmiBool QueryBalanceNow(int step);
  void work(CentralLB::LDStats* stats, int count);

};


#endif /* _ORBLB_H_ */

/*@}*/
