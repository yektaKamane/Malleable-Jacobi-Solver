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

#ifndef NBORBASELB_H
#define NBORBASELB_H

#include <BaseLB.h>
#include "NborBaseLB.decl.h"

void CreateNborBaseLB();

/// for backward compatibility
typedef LBMigrateMsg NLBMigrateMsg;

class NLBStatsMsg;

class NborBaseLB : public BaseLB
{
private:
  CProxy_NborBaseLB  thisProxy;
public:
  NborBaseLB(const CkLBOptions &);
  NborBaseLB(CkMigrateMessage *m):BaseLB(m) {}
  ~NborBaseLB();

  int useDefCtor(void){ return 1; }
  static void staticAtSync(void*);
  void AtSync(void); // Everything is at the PE barrier

  void ReceiveStats(NLBStatsMsg *); 		// Receive stats on PE 0
  void ResumeClients();
  void ReceiveMigration(LBMigrateMsg *); 	// Receive migration data

  // Migrated-element callback
  static void staticMigrated(void* me, LDObjHandle h);
  void Migrated(LDObjHandle h);

  void MigrationDone(void);  // Call when migration is complete
  int step() { return mystep; };

  struct LDStats {  // Passed to Strategy
    int from_pe;
    double total_walltime;
    double total_cputime;
    double idletime;
    double bg_walltime;
    double bg_cputime;
    double obj_walltime;
    double obj_cputime;
    int pe_speed;
    CmiBool available;
    CmiBool move;

    int n_objs;
    LDObjData* objData;
    int n_comm;
    LDCommData* commData;
  };

protected:
  virtual CmiBool QueryBalanceNow(int) { return CmiTrue; };  
  virtual CmiBool QueryMigrateStep(int) { return CmiTrue; };  
  virtual LBMigrateMsg* Strategy(LDStats* stats,int count);

  virtual int max_neighbors() {
    if (CkNumPes() > 2) return 2;
    else return (CkNumPes()-1);
  }

  virtual int num_neighbors() {
    if (CkNumPes() > 2) return 2;
    else return (CkNumPes()-1);
  };

  virtual void neighbors(int* _n) {
    _n[0] = (CkMyPe() + CkNumPes() -1) % CkNumPes();
    _n[1] = (CkMyPe() + 1) % CkNumPes();
  };

  int NeighborIndex(int pe);   // return the neighbor array index

  /*
  struct {
    int pe_speed;
    double total_walltime;
    double total_cputime;
    double idletime;
    double bg_walltime;
    double bg_cputime;
    int obj_data_sz;
    LDObjData* objData;
    int comm_data_sz;
    LDCommData* commData;
    double obj_walltime;
    double obj_cputime;
  } myStats;
  */
  LDStats myStats;

private:
  void FindNeighbors();
  NLBStatsMsg* AssembleStats();

  int mystep;
  int stats_msg_count;
  NLBStatsMsg** statsMsgsList;
  LDStats* statsDataList;
  int migrates_completed;
  int migrates_expected;
  LBMigrateMsg** mig_msgs;
  int mig_msgs_received;
  int mig_msgs_expected;
  int* neighbor_pes;
  int receive_stats_ready;
  double start_lb_time;
};

class NLBStatsMsg : public CMessage_NLBStatsMsg {
public:
  int from_pe;
  int serial;
  int pe_speed;
  double total_walltime;
  double total_cputime;
  double idletime;
  double bg_walltime;
  double bg_cputime;
  double obj_walltime;   // may not needed
  double obj_cputime;   // may not needed
  int n_objs;
  LDObjData *objData;
  int n_comm;
  LDCommData *commData;
}; 

#endif /* NBORBASELB_H */

/*@}*/
