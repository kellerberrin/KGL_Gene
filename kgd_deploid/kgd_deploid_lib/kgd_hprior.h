//
// Created by kellerberrin on 13/05/18.
//

#ifndef KGD_HPRIOR_H
#define KGD_HPRIOR_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_mersenne_twister.h"
#include "kgd_deploid_io.h"
#include "kgd_ibdconfig.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



class Hprior {

#ifdef UNITTEST
  friend class TestHprior;
  friend class TestMcmcMachinery;
  friend class TestIBDpath;
#endif

  friend class IBDpath;
  friend class McmcMachinery;
  friend class DEploidIO;

  Hprior() = default;
  ~Hprior() = default;


  IBDconfiguration ibdConfig;
  size_t kStrain_;
  size_t nLoci_;
  std::vector<double> plaf_;
  std::vector<std::vector<double> > priorProb; // size: nState x nLoci
  std::vector<std::vector<double> > priorProbTrans; // size: nLoci x nState
  std::vector<size_t> stateIdx; // size: nState
  std::vector<size_t> stateIdxFreq;

  std::vector<std::vector<int> > hSet; // size: nState x kStrain
  size_t nState_;

  void buildHprior(size_t kStrain, std::vector<double> &plaf);

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t getkStrain() const { return kStrain_; }

  void setnLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }

  void transposePriorProbs();

  size_t nState() const { return nState_; }

  std::vector<size_t> effectiveK;

  size_t nPattern() const { return effectiveK.size(); }

  std::vector<std::string> getIBDconfigureHeader();

};



}   // organization level namespace
}   // project level namespace


#endif //KGD_HPRIOR_H
