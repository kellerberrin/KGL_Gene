//
// Created by kellerberrin on 12/05/18.
//

#include "kgd_update_single_haplotype.h"
#include "kgl_exec_env.h"
#include <algorithm>    // std::reverse
#include <cstdlib>      // div


namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;


kgd::UpdateSingleHap::UpdateSingleHap(const std::vector<double> &refCount,
                                      const std::vector<double> &altCount,
                                      const std::vector<double> &plaf,
                                      const std::vector<double> &expectedWsaf,
                                      const std::vector<double> &proportion,
                                      const std::vector<std::vector<double> > &haplotypes,
                                      std::shared_ptr<RandomGenerator> randomGenerator,
                                      size_t segmentStartIndex,
                                      size_t nLoci,
                                      std::shared_ptr<Panel> panel,
                                      double missCopyProb,
                                      double scalingFactor,
                                      size_t strainIndex) : UpdateHap(refCount,
                                                                      altCount,
                                                                      expectedWsaf,
                                                                      plaf,
                                                                      proportion,
                                                                      haplotypes,
                                                                      randomGenerator,
                                                                      segmentStartIndex,
                                                                      nLoci,
                                                                      panel,
                                                                      missCopyProb,
                                                                      scalingFactor) {
  strainIndex_ = strainIndex;
  siteOfOneSwitchOne = std::vector<double>(nLoci);
  siteOfOneMissCopyOne = std::vector<double>(nLoci);

}


void kgd::UpdateSingleHap::core(const std::vector<double> &refCount,
                                const std::vector<double> &altCount,
                                const std::vector<double> &plaf,
                                const std::vector<double> &expectedWsaf,
                                const std::vector<double> &proportion,
                                const std::vector<std::vector<double> > &haplotypes) {

  calcExpectedWsaf(expectedWsaf, proportion, haplotypes);
  calcHapLLKs(refCount, altCount);

  if (panel_) {

    buildEmission(missCopyProb_);
    calcFwdProbs();
    samplePaths();
    addMissCopying(missCopyProb_);

  } else {

    sampleHapIndependently(plaf);

  }

  updateLLK();

}


void kgd::UpdateSingleHap::painting(std::vector<double> &refCount,
                                    std::vector<double> &altCount,
                                    std::vector<double> &expectedWsaf,
                                    std::vector<double> &proportion,
                                    std::vector<std::vector<double> > &haplotypes) {

  calcExpectedWsaf(expectedWsaf, proportion, haplotypes);
  calcHapLLKs(refCount, altCount);
  buildEmission(missCopyProb_);
  calcFwdBwdProbs();

}


void kgd::UpdateSingleHap::calcBwdProbs() {

  std::vector<double> bwdLast(nPanel_, 0.0);

  for (size_t i = 0; i < nPanel_; i++) {

    bwdLast[i] = 1.0;

  }

  Utility::normalizeBySum(bwdLast);
  assert(bwdProbs_.size() == 0);
  bwdProbs_.push_back(bwdLast);

  int j = (nLoci_ - 1);
  while (j > 0) {

    size_t hapIndexBack = segmentStartIndex_ + j;
    std::vector<double> bwdTmp(nPanel_, 1.0);
    double pRecEachHap = panel_->pRecEachHap_[hapIndexBack - 1];
    double pNoRec = panel_->pNoRec_[hapIndexBack - 1];

    for (size_t i = 0; i < nPanel_; i++) {

      bwdTmp[i] = 0.0;

      for (size_t ii = 0; ii < nPanel_; ii++) {

        bwdTmp[i] += emission_[j][panel_->content_[hapIndexBack][ii]] * bwdProbs_.back()[ii] * pRecEachHap;

        if (i == ii) {

          bwdTmp[i] += emission_[j][panel_->content_[hapIndexBack][ii]] * bwdProbs_.back()[ii] * pNoRec;

        }

      }

    }

    Utility::normalizeBySum(bwdTmp);
    bwdProbs_.push_back(bwdTmp);
    j--;

  }

  if (bwdProbs_.size() != nLoci_) {

    throw LociNumberUnequal("here");

  }

  assert (bwdProbs_.size() == nLoci_);

}

void kgd::UpdateSingleHap::calcFwdBwdProbs() {

  calcFwdProbs();
  calcBwdProbs();

  assert (fwdBwdProbs_.size() == 0);

  for (size_t j = 0; j < nLoci_; j++) {

    std::vector<double> fwdBwdTmp(nPanel_, 0.0);

    for (size_t i = 0; i < nPanel_; i++) {

      fwdBwdTmp[i] = fwdProbs_[j][i] * bwdProbs_[nLoci_ - j - 1][i];

    }

    Utility::normalizeBySum(fwdBwdTmp);
    fwdBwdProbs_.push_back(fwdBwdTmp);

  }

  assert (fwdBwdProbs_.size() == nLoci_);

}

void kgd::UpdateSingleHap::calcExpectedWsaf(const std::vector<double> &expectedWsaf,
                                            const std::vector<double> &proportion,
                                            const std::vector<std::vector<double> > &haplotypes) {

  assert (expectedWsaf0_.size() == 0);
  assert (expectedWsaf1_.size() == 0);

  expectedWsaf0_ = std::vector<double>(expectedWsaf.begin() + segmentStartIndex_,
                                       expectedWsaf.begin() + (segmentStartIndex_ + nLoci_));
  size_t hapIndex = segmentStartIndex_;

  for (size_t i = 0; i < expectedWsaf0_.size(); i++) {

    expectedWsaf0_[i] -= proportion[strainIndex_] * haplotypes[hapIndex][strainIndex_];

    assert (expectedWsaf0_[i] >= 0);

    assert (expectedWsaf0_[i] <= 1);
    hapIndex++;

  }

  expectedWsaf1_ = expectedWsaf0_;

  for (size_t i = 0; i < expectedWsaf1_.size(); i++) {

    expectedWsaf1_[i] += proportion[strainIndex_];
    assert (expectedWsaf1_[i] >= 0);

  }

  assert (expectedWsaf0_.size() == nLoci_);
  assert (expectedWsaf1_.size() == nLoci_);

}


void kgd::UpdateSingleHap::buildEmission(double missCopyProb) {

  std::vector<double> noMissProb(nLoci_, log(1.0 - missCopyProb));
  std::vector<double> t1omu = Utility::vecSum(llk0_, noMissProb); // t1 one minus u
  std::vector<double> t2omu = Utility::vecSum(llk1_, noMissProb); // t2 one minus u


  std::vector<double> missProb(nLoci_, log(missCopyProb));
  std::vector<double> t1u = Utility::vecSum(llk0_, missProb);
  std::vector<double> t2u = Utility::vecSum(llk1_, missProb);

  assert(emission_.size() == 0);

  for (size_t i = 0; i < nLoci_; i++) {

    std::vector<double> tmp({t1omu[i], t2omu[i], t1u[i], t2u[i]});
    double tmaxTmp = Utility::max_value(tmp);
    std::vector<double> emissRow({exp(t1omu[i] - tmaxTmp) + exp(t2u[i] - tmaxTmp),
                                  exp(t2omu[i] - tmaxTmp) + exp(t1u[i] - tmaxTmp)});

    emission_.push_back(emissRow);

  }

}


void kgd::UpdateSingleHap::buildEmissionBasicVersion(double missCopyProb) {

  assert(emission_.size() == 0);

  for (size_t i = 0; i < nLoci_; i++) {

    std::vector<double> emissRow({exp(llk0_[i]) * (1.0 - missCopyProb) + exp(llk1_[i]) * missCopyProb,
                                  exp(llk1_[i]) * (1.0 - missCopyProb) + exp(llk0_[i]) * missCopyProb});

    emission_.push_back(emissRow);

  }

}


void kgd::UpdateSingleHap::calcFwdProbs() {

  size_t hapIndex = segmentStartIndex_;

  assert (fwdProbs_.size() == 0);

  std::vector<double> fwd1st(nPanel_, 0.0);
  for (size_t i = 0; i < nPanel_; i++) {

    auto panel_index = panel_->content_[hapIndex][i];
    fwd1st[i] = emission_[0][panel_index];

  }

  Utility::normalizeBySum(fwd1st);
  fwdProbs_.push_back(fwd1st);

  for (size_t j = 1; j < nLoci_; j++) {

    double pRecEachHap = panel_->pRecEachHap_[hapIndex];
    double pNoRec = panel_->pNoRec_[hapIndex];

    hapIndex++;

    double massFromRec = Utility::sumOfVec(fwdProbs_.back()) * pRecEachHap;

    std::vector<double> fwdTmp(nPanel_, 0.0);

    for (size_t i = 0; i < nPanel_; i++) {

      fwdTmp[i] = emission_[j][panel_->content_[hapIndex][i]] * (fwdProbs_.back()[i] * pNoRec + massFromRec);

    }

    Utility::normalizeBySum(fwdTmp);

    fwdProbs_.push_back(fwdTmp);

  }

  assert (fwdProbs_.size() == nLoci_);

}


void kgd::UpdateSingleHap::calcHapLLKs(const std::vector<double> &refCount,
                                       const std::vector<double> &altCount) {

  llk0_ = Utility::calcLLKs(refCount, altCount, expectedWsaf0_, segmentStartIndex_, nLoci_, scalingFactor());
  llk1_ = Utility::calcLLKs(refCount, altCount, expectedWsaf1_, segmentStartIndex_, nLoci_, scalingFactor());

  assert(llk0_.size() == nLoci_);
  assert(llk1_.size() == nLoci_);

}


void kgd::UpdateSingleHap::samplePaths() {

  assert (path_.size() == 0);
  // Sample path at the last position

  size_t pathTmp = Utility::sampleIndexGivenProp(recombRg_, fwdProbs_.back());
  size_t contentIndex = segmentStartIndex_ + nLoci_ - 1;

  path_.push_back(panel_->content_[contentIndex][pathTmp]);

  for (size_t j = (nLoci_ - 1); j > 0; j--) {

    contentIndex--;
    double pRecEachHap = panel_->pRecEachHap_[contentIndex];
    double pNoRec = panel_->pNoRec_[contentIndex];

    size_t previous_site = j - 1;
    std::vector<double> previousDist = fwdProbs_[previous_site];

    std::vector<double> weightOfNoRecAndRec({previousDist[pathTmp] * pNoRec,
                                             Utility::sumOfVec(previousDist) * pRecEachHap});

    Utility::normalizeBySum(weightOfNoRecAndRec);

    if (Utility::sampleIndexGivenProp(recombRg_, weightOfNoRecAndRec) == (size_t) 1) { // Switch one

      pathTmp = Utility::sampleIndexGivenProp(recombLevel2Rg_, previousDist);
      siteOfOneSwitchOne[j] += 1.0;

    }

    path_.push_back(panel_->content_[contentIndex][pathTmp]);

  }

  reverse(path_.begin(), path_.end());
  assert(path_.size() == nLoci_);

}


void kgd::UpdateSingleHap::addMissCopying(double missCopyProb) {

  assert(hap_.size() == 0);

  for (size_t i = 0; i < nLoci_; i++) {

    double tmpMax = Utility::max_value(std::vector<double>({llk0_[i], llk1_[i]}));

    std::vector<double> emissionTmp({exp(llk0_[i] - tmpMax), exp(llk1_[i] - tmpMax)});

    std::vector<double> sameDiffDist({emissionTmp[path_[i]] * (1.0 - missCopyProb), // probability of the same
                                      emissionTmp[(size_t) (1 - path_[i])] * missCopyProb}); // probability of differ

    Utility::normalizeBySum(sameDiffDist);

    if (Utility::sampleIndexGivenProp(missCopyRg_, sameDiffDist) == 1) {

      hap_.push_back(1 - path_[i]); // differ
      siteOfOneMissCopyOne[i] += 1.0;

    } else {

      hap_.push_back(path_[i]); // same

    }

  }

  assert (hap_.size() == nLoci_);

}


void kgd::UpdateSingleHap::sampleHapIndependently(const std::vector<double> &plaf) {

  assert(hap_.size() == 0);

  size_t plafIndex = segmentStartIndex_;

  for (size_t i = 0; i < nLoci_; i++) {

    double tmpMax = Utility::max_value(std::vector<double>({llk0_[i], llk1_[i]}));
    std::vector<double> tmpDist({exp(llk0_[i] - tmpMax) * (1.0 - plaf[plafIndex]),
                                 exp(llk1_[i] - tmpMax) * plaf[plafIndex]});
    Utility::normalizeBySum(tmpDist);
    hap_.push_back((double) Utility::sampleIndexGivenProp(recombRg_, tmpDist));
    plafIndex++;

  }

  assert (hap_.size() == nLoci_);

}


void kgd::UpdateSingleHap::updateLLK() {

  newLLK = std::vector<double>(nLoci_, 0.0);

  for (size_t i = 0; i < nLoci_; i++) {
    if (hap_[i] == 0) {

      newLLK[i] = llk0_[i];

    } else if (hap_[i] == 1) {

      newLLK[i] = llk1_[i];

    } else {

      throw ShouldNotBeCalled();

    }

  }

}

