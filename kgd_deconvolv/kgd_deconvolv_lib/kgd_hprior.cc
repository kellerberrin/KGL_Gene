//
// Created by kellerberrin on 13/05/18.
//


#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "kgd_deconvolv_app.h"
#include "kgd_combinatorial.h"
#include "kgd_hprior.h"



namespace kgd = kellerberrin::deconvolv;


void kgd::Hprior::buildHprior(size_t kStrain, const std::vector<double> &plaf) {

  ibd_config_.buildIBDconfiguration(kStrain);

  effectiveK_ = ibd_config_.effectiveK();
  n_state_entries_ = 0;
  pop_allele_freq_ = plaf;
  setKstrain(kStrain);
  setnLoci(pop_allele_freq_.size());

  /// Get all the permutations (including 0, ..., 0) of k strains.
  BinaryPermutations binary_permutations(getkStrain());
  std::vector<std::vector<size_t> > hSetBase = binary_permutations.enumerateBinaryMatrix();

  size_t state_i = 0;

  std::vector<std::vector<double> > prior_probs; // size: (nState x entries) x nLoci; probabilities.

  /// For all states generated.
  for (auto state : ibd_config_.states()) {

    std::set<size_t> stateUnique(state.begin(), state.end());

    assert(stateUnique.size() == effectiveK_[state_i]);

    std::vector<std::vector<size_t> > hSetBaseTmp = hSetBase;

    for (size_t j = 0; j < getkStrain(); ++j) {

      for (size_t i = 0; i < hSetBase.size(); ++i) {

        hSetBaseTmp[i][j] = hSetBase[i][state[j]];

      }

    }

    std::vector<std::vector<size_t>> hSetBaseTmpUnique = Utility::uniqueMatrixColumns(hSetBaseTmp); // uu

//*****Debug.
    ExecEnv::log().info("hPrior::hSetBaseTmpUnique.size: {}", hSetBaseTmpUnique.size());
    for (auto colvec : hSetBaseTmpUnique) {

      std::stringstream ss;
      for (auto element : colvec) {

        ss << element << ", ";

      }

      ExecEnv::log().info("hPrior::hSetBaseTmpUnique: {} State: {}, EffectiveK: {}", ss.str(), state_i, effectiveK_[state_i]);

    }
//*****Debug.

    size_t sizeOfhSetBaseTmpUnique = hSetBaseTmpUnique.size();

    state_idx_freq_.push_back(sizeOfhSetBaseTmpUnique);

    for (size_t i = 0; i < sizeOfhSetBaseTmpUnique; i++) {

      size_t tmpSum = 0;

#define UNIQUE_PROB
#ifdef UNIQUE_PROB

      for (auto uniqSt : stateUnique) {

        tmpSum += hSetBaseTmpUnique[i][uniqSt];

      }

      size_t tmpDiff = stateUnique.size() - tmpSum;

#else

      for (auto value : hSetBaseTmpUnique[i]) {

        tmpSum += value;

      }

      size_t tmpDiff = hSetBaseTmpUnique[i].size() - tmpSum;

#endif

      std::vector<double> hPriorTmp;

      for (size_t site = 0; site < nLoci(); ++site) {

        double site_value = std::pow(pop_allele_freq_[site], static_cast<double>(tmpSum)) * std::pow((1.0 - pop_allele_freq_[site]), static_cast<double>(tmpDiff));

        hPriorTmp.push_back(site_value);

      }

//*****Debug.
      size_t prob_index = 250;

      ExecEnv::log().info("hPrior::hPriorTmp.size: {}, tmpSum: {}, tmpDiff: {}, state: {}, prob[{}]: {}, allele[{}]: {}, hSet_column.size :{}",
                          hPriorTmp.size(), tmpSum, tmpDiff, state_i, prob_index, hPriorTmp[prob_index], prob_index,
                          pop_allele_freq_[prob_index], hSetBaseTmpUnique[i].size());
//*****Debug.

      prior_probs.push_back(hPriorTmp);
      h_set_.push_back(hSetBaseTmpUnique[i]);

      n_state_entries_++;
      state_idx_.push_back(state_i);

    }

    state_i++;

  }

//******Debug
  ExecEnv::log().info("h_set_.size:{}, prior_probs.size: {}, prior_probs[0].size: {}, state_idx_.size: {}, state_idx_freq_.size: {},",
                      h_set_.size(), prior_probs.size(), prior_probs[0].size(), state_idx_.size(), state_idx_freq_.size());

//*******Debug

  generateEntryOffsets();

  transposePriorProbs(prior_probs);


}


void kgd::Hprior::generateEntryOffsets() {

  size_t state_entry_offset = 0;
  for (auto state_size : getStateIdxFreq()) {

    state_entry_offset_.push_back(state_entry_offset);

    state_entry_offset += state_size;

  }

}


void kgd::Hprior::transposePriorProbs(std::vector<std::vector<double> >& prior_probs) {

  prior_probs_trans_.clear();

  for (size_t site = 0; site < nLoci(); ++site) {

    std::vector<double> priorProbTransTmp(nStateEntries());

    for (size_t state_entry = 0; state_entry < nStateEntries(); ++state_entry) {

      priorProbTransTmp[state_entry] = prior_probs[state_entry][site];
    }

    prior_probs_trans_.push_back(priorProbTransTmp);

  }

}


std::vector<std::string> kgd::Hprior::getIBDconfigureHeader() const {

  return ibd_config_.getIBDconfigureHeader();

}

