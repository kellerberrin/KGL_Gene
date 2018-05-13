//
// Created by kellerberrin on 12/05/18.
//

#include "kgd_update_pair_haplotype.h"
#include "kgl_exec_env.h"
#include <algorithm>    // std::reverse
#include <cstdlib>      // div


namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;


kgd::UpdatePairHap::UpdatePairHap(std::vector<double> &refCount,
                                  std::vector<double> &altCount,
                                  std::vector<double> &plaf,
                                  std::vector<double> &expectedWsaf,
                                  std::vector<double> &proportion,
                                  std::vector<std::vector<double> > &haplotypes,
                                  std::shared_ptr<RandomGenerator> randomGenerator,
                                  size_t segmentStartIndex,
                                  size_t nLoci,
                                  std::shared_ptr<Panel> panel,
                                  double missCopyProb,
                                  double scalingFactor,
                                  bool forbidCopyFromSame,
                                  size_t strainIndex1,
                                  size_t strainIndex2) : UpdateHap(refCount,
                                                                   altCount,
                                                                   plaf,
                                                                   expectedWsaf,
                                                                   proportion,
                                                                   haplotypes,
                                                                   randomGenerator,
                                                                   segmentStartIndex,
                                                                   nLoci,
                                                                   panel,
                                                                   missCopyProb,
                                                                   scalingFactor) {
  strainIndex1_ = strainIndex1;
  strainIndex2_ = strainIndex2;
  forbidCopyFromSame_ = forbidCopyFromSame;
  siteOfTwoSwitchOne = std::vector<double>(nLoci);
  siteOfTwoMissCopyOne = std::vector<double>(nLoci);
  siteOfTwoSwitchTwo = std::vector<double>(nLoci);
  siteOfTwoMissCopyTwo = std::vector<double>(nLoci);

}


void kgd::UpdatePairHap::core(std::vector<double> &refCount,
                              std::vector<double> &altCount,
                              std::vector<double> &plaf,
                              std::vector<double> &expectedWsaf,
                              std::vector<double> &proportion,
                              std::vector<std::vector<double> > &haplotypes) {

  calcExpectedWsaf(expectedWsaf, proportion, haplotypes);
  calcHapLLKs(refCount, altCount);

  if (panel_ != NULL) {

    buildEmission(missCopyProb_);
    calcFwdProbs(forbidCopyFromSame_);
    samplePaths();
    addMissCopying(missCopyProb_);

  } else {

    sampleHapIndependently(plaf);

  }

  updateLLK();

}


void kgd::UpdatePairHap::calcExpectedWsaf(std::vector<double> &expectedWsaf,
                                          std::vector<double> &proportion,
                                          std::vector<std::vector<double> > &haplotypes) {


  expectedWsaf00_ = std::vector<double>(expectedWsaf.begin() + segmentStartIndex_,
                                        expectedWsaf.begin() + (segmentStartIndex_ + nLoci_));

  size_t hapIndex = segmentStartIndex_;

  for (size_t i = 0; i < expectedWsaf00_.size(); ++i) {

    double reduction = (proportion[strainIndex1_] * haplotypes[hapIndex][strainIndex1_]) + (proportion[strainIndex2_] * haplotypes[hapIndex][strainIndex2_]);

    double update = expectedWsaf00_[i] - reduction;

    if (update < 0 or update >= 1) {

      kgl::ExecEnv::log().error("file: {}, line: {}, bounds error not in [0, 1) update: {}, expectedWsaf00_[{}] = {}, reduction: {}",
                                __FILE__, __LINE__, update, i, expectedWsaf00_[i], reduction);

    }

    expectedWsaf00_[i] = update;

    assert (expectedWsaf00_[i] >= 0);
    assert (expectedWsaf00_[i] <= 1); // This was triggered by expectedWsaf00_[i] == 1.

    hapIndex++;

  }

  expectedWsaf10_ = expectedWsaf00_;
  for (size_t i = 0; i < expectedWsaf10_.size(); ++i) {

    expectedWsaf10_[i] += proportion[strainIndex1_];

  }

  expectedWsaf01_ = expectedWsaf00_;
  for (size_t i = 0; i < expectedWsaf01_.size(); ++i) {

    expectedWsaf01_[i] += proportion[strainIndex2_];

  }

  expectedWsaf11_ = expectedWsaf00_;
  for (size_t i = 0; i < expectedWsaf11_.size(); ++i) {

    expectedWsaf11_[i] += (proportion[strainIndex1_] + proportion[strainIndex2_]);

  }

}


void kgd::UpdatePairHap::calcHapLLKs(std::vector<double> &refCount, std::vector<double> &altCount) {

  llk00_ = calcLLKs(refCount, altCount, expectedWsaf00_, segmentStartIndex_, nLoci_, scalingFactor());
  llk10_ = calcLLKs(refCount, altCount, expectedWsaf10_, segmentStartIndex_, nLoci_, scalingFactor());
  llk01_ = calcLLKs(refCount, altCount, expectedWsaf01_, segmentStartIndex_, nLoci_, scalingFactor());
  llk11_ = calcLLKs(refCount, altCount, expectedWsaf11_, segmentStartIndex_, nLoci_, scalingFactor());
  assert(llk00_.size() == nLoci_);
  assert(llk10_.size() == nLoci_);
  assert(llk01_.size() == nLoci_);
  assert(llk11_.size() == nLoci_);

}


void kgd::UpdatePairHap::buildEmission(double missCopyProb) {

  std::vector<double> noMissProb(nLoci_, log(1.0 - missCopyProb));
  std::vector<double> missProb(nLoci_, log(missCopyProb));
  std::vector<double> noNo = vecSum(noMissProb, noMissProb);
  std::vector<double> misMis = vecSum(missProb, missProb);
  std::vector<double> misNo = vecSum(noMissProb, missProb);

  std::vector<double> tmp_00_1 = vecSum(llk00_, noNo);
  std::vector<double> tmp_00_2 = vecSum(llk10_, misNo);
  std::vector<double> tmp_00_3 = vecSum(llk01_, misNo);
  std::vector<double> tmp_00_4 = vecSum(llk11_, misMis);

  std::vector<double> tmp_01_1 = vecSum(llk01_, noNo);
  std::vector<double> tmp_01_2 = vecSum(llk00_, misNo);
  std::vector<double> tmp_01_3 = vecSum(llk11_, misNo);
  std::vector<double> tmp_01_4 = vecSum(llk10_, misMis);

  std::vector<double> tmp_10_1 = vecSum(llk10_, noNo);
  std::vector<double> tmp_10_2 = vecSum(llk00_, misNo);
  std::vector<double> tmp_10_3 = vecSum(llk11_, misNo);
  std::vector<double> tmp_10_4 = vecSum(llk01_, misMis);

  std::vector<double> tmp_11_1 = vecSum(llk11_, noNo);
  std::vector<double> tmp_11_2 = vecSum(llk10_, misNo);
  std::vector<double> tmp_11_3 = vecSum(llk01_, misNo);
  std::vector<double> tmp_11_4 = vecSum(llk00_, misMis);

  assert(emission_.size() == 0);

  for (size_t i = 0; i < nLoci_; i++) {

    std::vector<double> tmp({tmp_00_1[i], tmp_00_2[i], tmp_00_3[i], tmp_00_4[i],
                             tmp_01_1[i], tmp_01_2[i], tmp_01_3[i], tmp_01_4[i],
                             tmp_10_1[i], tmp_10_2[i], tmp_10_3[i], tmp_10_4[i],
                             tmp_11_1[i], tmp_11_2[i], tmp_11_3[i], tmp_11_4[i]});

    double tmaxTmp = max_value(tmp);

    std::vector<double> emissRow(
    {exp(tmp_00_1[i] - tmaxTmp) + exp(tmp_00_2[i] - tmaxTmp) + exp(tmp_00_3[i] - tmaxTmp) + exp(tmp_00_4[i] - tmaxTmp),
     exp(tmp_01_1[i] - tmaxTmp) + exp(tmp_01_2[i] - tmaxTmp) + exp(tmp_01_3[i] - tmaxTmp) + exp(tmp_01_4[i] - tmaxTmp),
     exp(tmp_10_1[i] - tmaxTmp) + exp(tmp_10_2[i] - tmaxTmp) + exp(tmp_10_3[i] - tmaxTmp) + exp(tmp_10_4[i] - tmaxTmp),
     exp(tmp_11_1[i] - tmaxTmp) + exp(tmp_11_2[i] - tmaxTmp) + exp(tmp_11_3[i] - tmaxTmp) +
     exp(tmp_11_4[i] - tmaxTmp)});

    emission_.push_back(emissRow);

  }

  assert(emission_.size() == nLoci_);

}


std::vector<double> kgd::UpdatePairHap::computeRowMarginalDist(std::vector<std::vector<double> > &probDist) { // Sum of Rows

  std::vector<double> marginalDist(probDist.size(), 0.0);

  for (size_t i = 0; i < probDist.size(); i++) {

    marginalDist[i] = sumOfVec(probDist[i]);

  }

  return marginalDist;

}


std::vector<double> kgd::UpdatePairHap::computeColMarginalDist(std::vector<std::vector<double> > &probDist) { // Sum of Cols

  std::vector<double> marginalDist(probDist.size(), 0.0);

  for (size_t coli = 0; coli < probDist[0].size(); coli++) {

    for (size_t rowi = 0; rowi < probDist.size(); rowi++) {

      marginalDist[coli] += probDist[rowi][coli];

    }

  }

  return marginalDist;

}


void kgd::UpdatePairHap::calcFwdProbs(bool forbidCopyFromSame) {

  size_t hapIndex = segmentStartIndex_;
  assert (fwdProbs_.size() == 0);
  std::vector<std::vector<double> > fwd1st;

  for (size_t i = 0; i < nPanel_; i++) { // Row of the matrix

    size_t rowObs = (size_t) panel_->content_[0][i];
    std::vector<double> fwd1stRow(nPanel_, 0.0);

    for (size_t ii = 0; ii < nPanel_; ii++) { // Column of the matrix

      if (forbidCopyFromSame && i == ii) continue;

      size_t colObs = (size_t) panel_->content_[hapIndex][ii];
      size_t obs = rowObs * 2 + colObs;
      fwd1stRow[ii] = emission_[0][obs];

    }

    fwd1st.push_back(fwd1stRow);

  }

  normalizeBySumMat(fwd1st);
  fwdProbs_.push_back(fwd1st);

  for (size_t j = 1; j < nLoci_; j++) {

    double recRec = panel_->pRecRec_[hapIndex];
    double recNorec = panel_->pRecNoRec_[hapIndex];
    double norecNorec = panel_->pNoRecNoRec_[hapIndex];
    hapIndex++;

    std::vector<double> marginalOfRows = computeRowMarginalDist(fwdProbs_.back());
    std::vector<double> marginalOfCols = computeColMarginalDist(fwdProbs_.back());

    std::vector<std::vector<double> > fwdTmp;
    for (size_t i = 0; i < nPanel_; i++) {

      size_t rowObs = (size_t) panel_->content_[hapIndex][i];
      std::vector<double> fwdTmpRow(nPanel_, 0.0);
      for (size_t ii = 0; ii < nPanel_; ii++) {
        if (forbidCopyFromSame && i == ii) continue;

        size_t colObs = (size_t) panel_->content_[hapIndex][ii];
        size_t obs = rowObs * 2 + colObs;
        fwdTmpRow[ii] = emission_[j][obs] * (sumOfMat(fwdProbs_.back()) * recRec +
                                             fwdProbs_.back()[i][ii] * norecNorec +
                                             recNorec * (marginalOfRows[ii] + marginalOfCols[i]));

      }

      fwdTmp.push_back(fwdTmpRow);

    }

    normalizeBySumMat(fwdTmp);
    fwdProbs_.push_back(fwdTmp);

  }

}


std::vector<size_t> kgd::UpdatePairHap::sampleMatrixIndex(std::vector<std::vector<double> > &probDist) {

  size_t tmp = sampleIndexGivenProp(recombLevel2Rg_, reshapeMatToVec(probDist));
  div_t divresult;
  divresult = div((int) tmp, (int) nPanel_);

  return std::vector<size_t>({(size_t) divresult.quot, (size_t) divresult.rem});

}


void kgd::UpdatePairHap::samplePaths() {

  assert (path1_.size() == 0);
  assert (path2_.size() == 0);

  std::vector<size_t> tmpPath = sampleMatrixIndex(fwdProbs_[nLoci_ - 1]);
  size_t rowI = tmpPath[0];
  size_t colJ = tmpPath[1];
  size_t contentIndex = segmentStartIndex_ + nLoci_ - 1;

  path1_.push_back(panel_->content_[contentIndex][rowI]);
  path2_.push_back(panel_->content_[contentIndex][colJ]);

  for (size_t j = (nLoci_ - 1); j > 0; j--) {

    --contentIndex;
    double recRec = panel_->pRecRec_[contentIndex];
    double recNorec = panel_->pRecNoRec_[contentIndex];
    double norecNorec = panel_->pNoRecNoRec_[contentIndex];

    size_t previous_site = j - 1;
    std::vector<std::vector<double> > previousDist = fwdProbs_[previous_site];
    double previousProbij = previousDist[rowI][colJ];

    std::vector<double> rowIdist = previousDist[rowI];
    double tmpRowSum = sumOfVec(rowIdist);

    std::vector<double> colJdist;
    for (auto const &array: previousDist) {

      colJdist.push_back(array[colJ]);

    }

    assert(nPanel_ == colJdist.size());

    double tmpColSum = sumOfVec(colJdist);

    std::vector<double> weightOfFourCases(
    {recRec * sumOfMat(previousDist),           // recombination happened on both strains
     recNorec * tmpRowSum,  // first strain no recombine, second strain recombine
     recNorec * tmpColSum,  // first strain recombine, second strain no recombine
     norecNorec * previousProbij}); // no recombine on either strain

    normalizeBySum(weightOfFourCases);

    size_t tmpCase = sampleIndexGivenProp(recombRg_, weightOfFourCases);

    if (tmpCase == (size_t) 0) { // switching both strains

      siteOfTwoSwitchTwo[j] += 1.0;
      tmpPath = sampleMatrixIndex(previousDist);
      rowI = tmpPath[0];
      colJ = tmpPath[1];
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else if (tmpCase == 1) { // switching second strain

      siteOfTwoSwitchOne[j] += 0.5;
      normalizeBySum(rowIdist);
      colJ = sampleIndexGivenProp(recombLevel2Rg_, rowIdist);
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else if (tmpCase == (size_t) 2) { // switching first strain

      siteOfTwoSwitchOne[j] += 0.5;
      normalizeBySum(colJdist);
      rowI = sampleIndexGivenProp(recombLevel2Rg_, colJdist);
      colJ = colJ;
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else if (tmpCase == 3) { // no switching

      rowI = rowI;
      colJ = colJ;
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else {

      throw ShouldNotBeCalled();

    }

    path1_.push_back(panel_->content_[contentIndex][rowI]);
    path2_.push_back(panel_->content_[contentIndex][colJ]);

  }

  reverse(path1_.begin(), path1_.end());
  reverse(path2_.begin(), path2_.end());

  assert(path1_.size() == nLoci_);
  assert(path2_.size() == nLoci_);

}


void kgd::UpdatePairHap::addMissCopying(double missCopyProb) {

  assert(hap1_.size() == 0);
  assert(hap2_.size() == 0);

  for (size_t i = 0; i < nLoci_; i++) {

    double tmpMax = max_value(std::vector<double>({llk00_[i], llk01_[i], llk10_[i], llk11_[i]}));

    std::vector<double> emissionTmp(
    {exp(llk00_[i] - tmpMax), exp(llk01_[i] - tmpMax), exp(llk10_[i] - tmpMax), exp(llk11_[i] - tmpMax)});

    std::vector<double> casesDist({emissionTmp[(size_t) (2 * path1_[i] + path2_[i])] * (1.0 - missCopyProb) *
                                   (1.0 - missCopyProb), // probability of both same
                                   emissionTmp[(size_t) (2 * path1_[i] + (1 - path2_[i]))] * (1.0 - missCopyProb) *
                                   missCopyProb,         // probability of same1diff2
                                   emissionTmp[(size_t) (2 * (1 - path1_[i]) + path2_[i])] * missCopyProb *
                                   (1.0 - missCopyProb),         // probability of same2diff1
                                   emissionTmp[(size_t) (2 * (1 - path1_[i]) + (1 - path2_[i]))] * missCopyProb *
                                   missCopyProb});              // probability of both differ
    normalizeBySum(casesDist);

    size_t tmpCase = sampleIndexGivenProp(missCopyRg_, casesDist);

    if (tmpCase == 0) {

      hap1_.push_back(path1_[i]);
      hap2_.push_back(path2_[i]);

    } else if (tmpCase == 1) {

      siteOfTwoMissCopyOne[i] += 0.5;
      hap1_.push_back(path1_[i]);
      hap2_.push_back(1.0 - path2_[i]);

    } else if (tmpCase == 2) {

      siteOfTwoMissCopyOne[i] += 0.5;
      hap1_.push_back(1.0 - path1_[i]);
      hap2_.push_back(path2_[i]);

    } else if (tmpCase == 3) {

      siteOfTwoMissCopyTwo[i] += 1.0;
      hap1_.push_back(1.0 - path1_[i]);
      hap2_.push_back(1.0 - path2_[i]);

    } else {

      throw ShouldNotBeCalled();

    }

  }

  assert (hap1_.size() == nLoci_);
  assert (hap2_.size() == nLoci_);

}


void kgd::UpdatePairHap::sampleHapIndependently(std::vector<double> &plaf) {

  assert(hap1_.size() == 0);
  assert(hap2_.size() == 0);

  size_t plafIndex = segmentStartIndex_;

  for (size_t i = 0; i < nLoci_; i++) {

    double tmpMax = max_value(std::vector<double>({llk00_[i], llk01_[i], llk10_[i], llk11_[i]}));

    std::vector<double> tmpDist({exp(llk00_[i] - tmpMax) * (1.0 - plaf[plafIndex]) * (1.0 - plaf[plafIndex]),
                                 exp(llk01_[i] - tmpMax) * (1.0 - plaf[plafIndex]) * plaf[plafIndex],
                                 exp(llk10_[i] - tmpMax) * (1.0 - plaf[plafIndex]) * plaf[plafIndex],
                                 exp(llk11_[i] - tmpMax) * plaf[plafIndex] * plaf[plafIndex]});

    normalizeBySum(tmpDist);

    size_t tmpCase = sampleIndexGivenProp(recombRg_, tmpDist);

    if (tmpCase == 0) {

      hap1_.push_back(0.0);
      hap2_.push_back(0.0);

    } else if (tmpCase == 1) {

      hap1_.push_back(0.0);
      hap2_.push_back(1.0);

    } else if (tmpCase == 2) {

      hap1_.push_back(1.0);
      hap2_.push_back(0.0);

    } else if (tmpCase == 3) {

      hap1_.push_back(1.0);
      hap2_.push_back(1.0);

    } else {

      throw ShouldNotBeCalled();

    }

    plafIndex++;

  }

  assert (hap1_.size() == nLoci_);
  assert (hap2_.size() == nLoci_);

}


void kgd::UpdatePairHap::updateLLK() {

  newLLK = std::vector<double>(nLoci_, 0.0);

  for (size_t i = 0; i < nLoci_; i++) {

    if (hap1_[i] == 0 && hap2_[i] == 0) {

      newLLK[i] = llk00_[i];

    } else if (hap1_[i] == 0 && hap2_[i] == 1) {

      newLLK[i] = llk01_[i];

    } else if (hap1_[i] == 1 && hap2_[i] == 0) {

      newLLK[i] = llk10_[i];

    } else if (hap1_[i] == 1 && hap2_[i] == 1) {

      newLLK[i] = llk11_[i];

    } else {

      throw ShouldNotBeCalled();

    }

  }

}
