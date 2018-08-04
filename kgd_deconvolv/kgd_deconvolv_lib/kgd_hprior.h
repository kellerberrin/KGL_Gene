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
namespace deconvolv {          // project level namespace



class Hprior {


public:

  Hprior() = default;
  ~Hprior() = default;

  const std::vector<std::vector<size_t> >& gethSet() const { return h_set_; }
  const std::vector<size_t>& getStateIdx() const { return state_idx_; }
  size_t nStates() const { return effectiveK_.size(); }
  size_t nStateEntries() const { return n_state_entries_; }
  const std::vector<size_t>& getStateEntryOffsets() const { return state_entry_offset_; }
  const std::vector<size_t>& getStateIdxFreq() const { return state_idx_freq_; }
  const std::vector<std::vector<double> >& getPriorProb() const { return prior_probs_trans_; }
  const std::vector<size_t>& getEffectiveK() const { return effectiveK_; }


  std::vector<std::string> getIBDconfigureHeader() const;

  void initializeHprior(size_t kStrain, const std::vector<double> &plaf) { buildHprior(kStrain, plaf); }

private:

  IBDconfiguration ibd_config_;
  size_t kStrain_;
  size_t nLoci_;
  std::vector<double> pop_allele_freq_; // nLoci
  std::vector<std::vector<double> > prior_probs_trans_; // size: nLoci x (nState x entries)
  std::vector<size_t> state_idx_; // size: (nState x entries) // state index.
  std::vector<size_t> state_idx_freq_;  // size: nstate (2^(effectiveK[state])
  std::vector<size_t> state_entry_offset_;  // size nStates, the offset of stateEntries for each state.
  std::vector<std::vector<size_t> > h_set_; // size: (nState x entries) x kStrain
  size_t n_state_entries_; // (nState x entries) increment from 0
  std::vector<size_t> effectiveK_; // nState

  void buildHprior(size_t kStrain, const std::vector<double> &plaf);

  void transposePriorProbs(std::vector<std::vector<double> >& prior_probs);

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t getkStrain() const { return kStrain_; }

  void setnLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }

  void generateEntryOffsets();

  };



}   // organization level namespace
}   // project level namespace


#endif //KGD_HPRIOR_H
