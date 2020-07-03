#ifndef KGROUP_H
#define KGROUP_H

#include "MGroup.h"
#include "ECMGroup.h"

/// An AZURE \f$ s,s' \f$ group.

/*!
 * In R-Matrix formalism, the equations required to calculate the cross section usually nested 
 * inside sums over entrance and exit channel spins.  For this reason AZURE groups reaction pathways 
 * according to their entrance and exit channel spins.  Each KGroup object is a container for vectors
 * of MGroup and ECMGroup objects.  
 */

class KGroup {
 public:
  KGroup(double, double);
  int NumMGroups() const;
  int NumECMGroups() const;
  int IsMGroup(MGroup);
  double GetS() const;
  double GetSp() const;
  void AddMGroup(MGroup);
  void AddECMGroup(ECMGroup);
  MGroup *GetMGroup(int);
  ECMGroup *GetECMGroup(int);
 private:
  double s_;
  double sp_;
  std::vector<MGroup> mgroups_;
  std::vector<ECMGroup> ec_mgroups_;
};

#endif
