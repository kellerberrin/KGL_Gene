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

public:

  Hprior() = default;
  ~Hprior() = default;

  const std::vector<std::vector<int> >& gethSet() const { return h_set_; }
  const std::vector<size_t>& getStateIdx() const { return state_idx_; }
  size_t nPattern() const { return effectiveK.size(); }
  size_t nState() const { return nState_; }
  const std::vector<size_t>& getStateIdxFreq() const { return state_idx_freq_; }
  const std::vector<std::vector<double> >& getPriorProb() const { return prior_probs_; }
  const std::vector<std::vector<double> >& getPriorProbTrans() const { return prior_probs_trans_; }
  const std::vector<size_t>& getEffectiveK() const { return effectiveK; }

  std::vector<std::string> getIBDconfigureHeader() const;

  void initializeHprior(size_t kStrain, const std::vector<double> &plaf) { buildHprior(kStrain, plaf); transposePriorProbs(); }

private:

  IBDconfiguration ibd_config_;
  size_t kStrain_;
  size_t nLoci_;
  std::vector<double> pop_allele_freq_;
  std::vector<std::vector<double> > prior_probs_; // size: nState x nLoci
  std::vector<std::vector<double> > prior_probs_trans_; // size: nLoci x nState
  std::vector<size_t> state_idx_; // size: nState
  std::vector<size_t> state_idx_freq_;
  std::vector<std::vector<int> > h_set_; // size: nState x kStrain
  size_t nState_;
  std::vector<size_t> effectiveK;

  void buildHprior(size_t kStrain, const std::vector<double> &plaf);

  void transposePriorProbs();

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t getkStrain() const { return kStrain_; }

  void setnLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }


};



}   // organization level namespace
}   // project level namespace


#endif //KGD_HPRIOR_H
