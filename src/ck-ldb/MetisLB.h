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

#ifndef _METISLB_H_
#define _METISLB_H_

#include "CentralLB.h"
#include "MetisLB.decl.h"

void CreateMetisLB();

class MetisLB : public CentralLB {
public:
  MetisLB(const CkLBOptions &);
  MetisLB(CkMigrateMessage *m):CentralLB(m) {}
private:
  CmiBool QueryBalanceNow(int step);
  void work(CentralLB::LDStats* stats, int count);
};

#define WEIGHTED 1
#define MULTI_CONSTRAINT 2

#endif /* _METISLB_H_ */

/*@}*/
