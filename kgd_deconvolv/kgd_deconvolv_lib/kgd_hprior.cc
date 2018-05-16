//
// Created by kellerberrin on 13/05/18.
//


#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "kgd_deconvolv_app.h"
#include "kgd_hprior.h"



namespace kgd = kellerberrin::deconvolv;


void kgd::Hprior::buildHprior(size_t kStrain, const std::vector<double> &plaf) {

  ibd_config_.buildIBDconfiguration(kStrain);

  effectiveK = ibd_config_.effectiveK();
  nState_ = 0;
  pop_allele_freq_ = plaf;
  setKstrain(kStrain);
  setnLoci(pop_allele_freq_.size());

  /// Get all the permutations (including 0, ..., 0) of k strains.
  std::vector<std::vector<int> > hSetBase = IBDconfiguration::enumerateBinaryMatrixOfK(getkStrain());
  size_t stateI = 0;

  /// For all states generated.
  for (auto state : ibd_config_.states()) {

    std::set<int> stateUnique(state.begin(), state.end());

    assert(stateUnique.size() == effectiveK[stateI]);

    std::vector<std::vector<int> > hSetBaseTmp = hSetBase;

    for (size_t j = 0; j < (size_t) getkStrain(); j++) {

      for (size_t i = 0; i < hSetBase.size(); i++) {

        hSetBaseTmp[i][j] = hSetBase[i][state[j]];

      }

    }

    std::vector<std::vector<int> > hSetBaseTmpUnique = IBDconfiguration::unique(hSetBaseTmp); // uu

    size_t sizeOfhSetBaseTmpUnique = hSetBaseTmpUnique.size();

    state_idx_freq_.push_back(sizeOfhSetBaseTmpUnique);
    //h.prior.i<-array(0, c(size.h.set.i, n.loci));

    for (size_t i = 0; i < sizeOfhSetBaseTmpUnique; i++) {
      //vector<int> hSetBaseTmpUniqueSubSet(); // uu[i,a.u,drop=F]
      int tmpSum = 0;

      for (int uniqSt : stateUnique) {

        tmpSum += hSetBaseTmpUnique[i][uniqSt];

      }
      // sumOfVec(hSetBaseTmpUnique[i]);
      int tmpDiff = stateUnique.size() - tmpSum;
      //cout << stateUnique.size() << " " << tmpSum << " " << tmpDiff<< endl;
      std::vector<double> hPriorTmp(nLoci());

      for (size_t site = 0; site < nLoci(); site++) {
        //cout << (1.0-pop_allele_freq_[site]) << " " << tmpDiff << " " <<pow((1.0-pop_allele_freq_[site]),(double)tmpDiff) << endl;

        hPriorTmp[site] = pow(pop_allele_freq_[site], (double) tmpSum) * pow((1.0 - pop_allele_freq_[site]), (double) tmpDiff);

      }

      prior_probs_.push_back(hPriorTmp);
      h_set_.push_back(hSetBaseTmpUnique[i]);

      nState_++;
      state_idx_.push_back(stateI);

    }

    stateI++;

  }

}


void kgd::Hprior::transposePriorProbs() {

  assert(prior_probs_trans_.size() == 0);

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> priorProbTransTmp(nState());

    for (size_t j = 0; j < nState(); j++) {

      priorProbTransTmp[j] = prior_probs_[j][i];
      //cout << prior_probs_[j][i] << " ";
    }

    prior_probs_trans_.push_back(priorProbTransTmp);
    //cout << endl;
  }

}


std::vector<std::string> kgd::Hprior::getIBDconfigureHeader() const {

  return ibd_config_.getIBDconfigureHeader();

}

