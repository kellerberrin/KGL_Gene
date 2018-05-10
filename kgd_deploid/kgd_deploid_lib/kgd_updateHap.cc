/*
 * kgd_deploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deploid.
 *
 * kgd_deploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "kgl_exec_env.h"
#include "kgd_updateHap.h"
#include <algorithm>    // std::reverse
#include <cstdlib>      // div

namespace kgl = kellerberrin::genome;


UpdateHap::UpdateHap( std::vector <double> &refCount,
                      std::vector <double> &altCount,
                      std::vector <double> &plaf,
                      std::vector <double> &expectedWsaf,
                      std::vector <double> &proportion,
                      std::vector < std::vector <double> > &haplotypes,
                      RandomGenerator* rg,
                      size_t segmentStartIndex,
                      size_t nLoci,
                      std::shared_ptr<Panel> panel,
                      double missCopyProb,
                      double scalingFactor) {

  panel_ = panel;
  nPanel_ = 0; // Initialize when panel is not given

  if (panel_) {

    setPanelSize( panel_->truePanelSize() );

  }

  kStrain_ = proportion.size();
  missCopyProb_ = missCopyProb;
  setScalingFactor(scalingFactor);
  recombRg_ = rg;
  recombLevel2Rg_ = rg;
  missCopyRg_ = rg;
  segmentStartIndex_ = segmentStartIndex;
  nLoci_ = nLoci;

}

void UpdateHap::core(std::vector <double> &refCount,
                           std::vector <double> &altCount,
                           std::vector <double> &plaf,
                           std::vector <double> &expectedWsaf,
                           std::vector <double> &proportion,
                           std::vector < std::vector <double> > &haplotypes ){ throw VirtualFunctionShouldNotBeCalled();};
void UpdateHap::calcExpectedWsaf( std::vector <double> & expectedWsaf, std::vector <double> &proportion, std::vector < std::vector <double> > &haplotypes){ throw VirtualFunctionShouldNotBeCalled(); };
void UpdateHap::calcHapLLKs( std::vector <double> &refCount, std::vector <double> &altCount){ throw VirtualFunctionShouldNotBeCalled();};
void UpdateHap::buildEmission( double missCopyProb ){ throw VirtualFunctionShouldNotBeCalled();};
void UpdateHap::samplePaths(){ throw VirtualFunctionShouldNotBeCalled();};
void UpdateHap::addMissCopying( double missCopyProb ){ throw VirtualFunctionShouldNotBeCalled();};
void UpdateHap::updateLLK(){ throw VirtualFunctionShouldNotBeCalled();};
void UpdateHap::sampleHapIndependently(std::vector <double> &plaf){ throw VirtualFunctionShouldNotBeCalled();};



UpdateSingleHap::UpdateSingleHap( std::vector <double> &refCount,
                                  std::vector <double> &altCount,
                                  std::vector <double> &plaf,
                                  std::vector <double> &expectedWsaf,
                                  std::vector <double> &proportion,
                                  std::vector < std::vector <double> > &haplotypes,
                                  RandomGenerator* rg,
                                  size_t segmentStartIndex,
                                  size_t nLoci,
                                  std::shared_ptr<Panel> panel, double missCopyProb, double scalingFactor,
                                  size_t strainIndex ): UpdateHap(refCount,
                                                                  altCount,
                                                                  expectedWsaf,
                                                                  plaf,
                                                                  proportion,
                                                                  haplotypes,
                                                                  rg,
                                                                  segmentStartIndex,
                                                                  nLoci,
                                                                  panel,
                                                                  missCopyProb,
                                                                  scalingFactor) {
  strainIndex_ = strainIndex;
  siteOfOneSwitchOne = std::vector <double>(nLoci);
  siteOfOneMissCopyOne = std::vector <double>(nLoci);

}


void UpdateSingleHap::core(std::vector <double> &refCount,
                           std::vector <double> &altCount,
                           std::vector <double> &plaf,
                           std::vector <double> &expectedWsaf,
                           std::vector <double> &proportion,
                           std::vector < std::vector <double> > &haplotypes ){

  calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
  calcHapLLKs(refCount, altCount);

  if (panel_){

    buildEmission( missCopyProb_ );
    calcFwdProbs();
    samplePaths();
    addMissCopying( missCopyProb_ );

  } else {

    sampleHapIndependently(plaf);

  }

  updateLLK();

}


void UpdateSingleHap::painting( std::vector <double> &refCount,
                                std::vector <double> &altCount,
                                std::vector <double> &expectedWsaf,
                                std::vector <double> &proportion,
                                std::vector < std::vector <double> > &haplotypes ){

  calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
  calcHapLLKs(refCount, altCount);
  buildEmission( missCopyProb_ );
  calcFwdBwdProbs();

}


void UpdateSingleHap::calcBwdProbs(){

  std::vector <double> bwdLast (nPanel_, 0.0);

  for ( size_t i = 0 ; i < nPanel_; i++){

    bwdLast[i] = 1.0;

  }

  normalizeBySum(bwdLast);
  assert(bwdProbs_.size() == 0 );
  bwdProbs_.push_back(bwdLast);

  int j = (nLoci_- 1);
  while (j > 0 ){

    size_t hapIndexBack = segmentStartIndex_ + j;
    std::vector <double> bwdTmp (nPanel_, 1.0);
    double pRecEachHap = panel_->pRecEachHap_[hapIndexBack-1];
    double pNoRec = panel_->pNoRec_[hapIndexBack-1];

    for ( size_t i = 0 ; i < nPanel_; i++) {

      bwdTmp[i] = 0.0;

      for ( size_t ii = 0 ; ii < nPanel_; ii++) {

        bwdTmp[i] += emission_[j][panel_->content_[hapIndexBack][ii]] * bwdProbs_.back()[ii] * pRecEachHap;

        if ( i == ii) {

          bwdTmp[i] += emission_[j][panel_->content_[hapIndexBack][ii]] * bwdProbs_.back()[ii] * pNoRec;

        }

      }

    }

    normalizeBySum(bwdTmp);
    bwdProbs_.push_back(bwdTmp);
    j--;

  }

  if (bwdProbs_.size() != nLoci_){

    throw LociNumberUnequal("here");

  }

  assert ( bwdProbs_.size() == nLoci_ );

}

void UpdateSingleHap::calcFwdBwdProbs(){

  calcFwdProbs();
  calcBwdProbs();

  assert (fwdBwdProbs_.size() == 0);

  for ( size_t j = 0; j < nLoci_; j++ ) {

    std::vector <double> fwdBwdTmp (nPanel_, 0.0);

    for ( size_t i = 0 ; i < nPanel_; i++ ){

      fwdBwdTmp[i] = fwdProbs_[j][i] * bwdProbs_[nLoci_-j-1][i];

    }

    normalizeBySum(fwdBwdTmp);
    fwdBwdProbs_.push_back(fwdBwdTmp);

  }

  assert (fwdBwdProbs_.size() == nLoci_ );

}

void UpdateSingleHap::calcExpectedWsaf( std::vector <double> & expectedWsaf, std::vector <double> &proportion, std::vector < std::vector <double> > &haplotypes ){

  assert ( expectedWsaf0_.size() == 0);
  assert ( expectedWsaf1_.size() == 0);

  expectedWsaf0_ = std::vector <double> (expectedWsaf.begin()+segmentStartIndex_, expectedWsaf.begin()+(segmentStartIndex_+nLoci_));
  size_t hapIndex = segmentStartIndex_;

  for ( size_t i = 0; i < expectedWsaf0_.size(); i++ ) {

    expectedWsaf0_[i] -= proportion[strainIndex_] * haplotypes[hapIndex][strainIndex_];

    assert (expectedWsaf0_[i] >= 0 );

    assert (expectedWsaf0_[i] <= 1 );
    hapIndex++;

  }

  expectedWsaf1_ = expectedWsaf0_;

  for ( size_t i = 0; i < expectedWsaf1_.size(); i++ ){

    expectedWsaf1_[i] += proportion[strainIndex_] ;
    assert (expectedWsaf1_[i] >= 0 );

  }

  assert ( expectedWsaf0_.size() == nLoci_ );
  assert ( expectedWsaf1_.size() == nLoci_ );

}


void UpdateSingleHap::buildEmission( double missCopyProb ) {

  std::vector <double> noMissProb (nLoci_, log(1.0 - missCopyProb));
  std::vector <double> t1omu = vecSum(llk0_, noMissProb); // t1 one minus u
  std::vector <double> t2omu = vecSum(llk1_, noMissProb); // t2 one minus u


  std::vector <double> missProb (nLoci_, log(missCopyProb));
  std::vector <double> t1u = vecSum(llk0_, missProb);
  std::vector <double> t2u = vecSum(llk1_, missProb);

  assert(emission_.size() == 0 );

  for ( size_t i = 0; i < nLoci_; i++){

    std::vector <double> tmp ({t1omu[i], t2omu[i], t1u[i], t2u[i]});
    double tmaxTmp = max_value(tmp);
    std::vector <double> emissRow ({exp(t1omu[i] - tmaxTmp) + exp(t2u[i] - tmaxTmp),
                               exp(t2omu[i] - tmaxTmp) + exp(t1u[i] - tmaxTmp)});

    emission_.push_back(emissRow);

  }

}


void UpdateSingleHap::buildEmissionBasicVersion( double missCopyProb ) {

  assert(emission_.size() == 0 );

  for ( size_t i = 0; i < nLoci_; i++) {

    std::vector <double> emissRow ({exp(llk0_[i])*(1.0-missCopyProb) + exp(llk1_[i])*missCopyProb,
                               exp(llk1_[i])*(1.0-missCopyProb) + exp(llk0_[i])*missCopyProb});

    emission_.push_back(emissRow);

  }

}


void UpdateSingleHap::calcFwdProbs(){

  size_t hapIndex = segmentStartIndex_;

  assert ( fwdProbs_.size() == 0 );

  std::vector <double> fwd1st (nPanel_, 0.0);
  for ( size_t i = 0 ; i < nPanel_; i++){

    auto panel_index = panel_->content_[hapIndex][i];
    fwd1st[i] = emission_[0][panel_index];

  }

  normalizeBySum(fwd1st);
  fwdProbs_.push_back(fwd1st);

  for ( size_t j = 1; j < nLoci_; j++ ){

    double pRecEachHap = panel_->pRecEachHap_[hapIndex];
    double pNoRec = panel_->pNoRec_[hapIndex];

    hapIndex++;

    double massFromRec = sumOfVec(fwdProbs_.back()) * pRecEachHap;

    std::vector <double> fwdTmp (nPanel_, 0.0);

    for ( size_t i = 0 ; i < nPanel_; i++){

        fwdTmp[i] = emission_[j][panel_->content_[hapIndex][i]] * (fwdProbs_.back()[i] * pNoRec + massFromRec);

    }

    normalizeBySum(fwdTmp);

    fwdProbs_.push_back(fwdTmp);

  }

  assert ( fwdProbs_.size() == nLoci_ );

}


void UpdateSingleHap::calcHapLLKs( std::vector <double> &refCount,
                                   std::vector <double> &altCount) {

  llk0_ = calcLLKs( refCount, altCount, expectedWsaf0_, segmentStartIndex_, nLoci_, scalingFactor() );
  llk1_ = calcLLKs( refCount, altCount, expectedWsaf1_, segmentStartIndex_, nLoci_, scalingFactor() );

  assert( llk0_.size() == nLoci_ );
  assert( llk1_.size() == nLoci_ );

}


void UpdateSingleHap::samplePaths(){

  assert ( path_.size() == 0 );
  // Sample path at the last position

  size_t pathTmp = sampleIndexGivenProp ( recombRg_, fwdProbs_.back() );
  size_t contentIndex = segmentStartIndex_ + nLoci_ - 1;

  path_.push_back( panel_->content_[contentIndex][pathTmp]);

  for ( size_t j = (nLoci_ - 1) ; j > 0; j-- ) {

    contentIndex--;
    double pRecEachHap = panel_->pRecEachHap_[contentIndex];
    double pNoRec = panel_->pNoRec_[contentIndex];

    size_t previous_site = j - 1;
    std::vector <double> previousDist = fwdProbs_[previous_site];

    std::vector <double> weightOfNoRecAndRec ({ previousDist[pathTmp]*pNoRec,
                                           sumOfVec(previousDist)*pRecEachHap});

    normalizeBySum(weightOfNoRecAndRec);

    if ( sampleIndexGivenProp(recombRg_, weightOfNoRecAndRec) == (size_t)1 ){ // Switch one

      pathTmp = sampleIndexGivenProp( recombLevel2Rg_, previousDist );
      siteOfOneSwitchOne[j] += 1.0;

    }

    path_.push_back(panel_->content_[contentIndex][pathTmp]);

  }

  reverse(path_.begin(), path_.end());
  assert(path_.size() == nLoci_);

}


void UpdateSingleHap::addMissCopying( double missCopyProb ){

  assert( hap_.size() == 0 );

  for ( size_t i = 0; i < nLoci_; i++) {

    double tmpMax = max_value ( std::vector <double>({llk0_[i], llk1_[i]}));

    std::vector <double> emissionTmp ({exp(llk0_[i]-tmpMax), exp(llk1_[i]-tmpMax)});

    std::vector <double> sameDiffDist ({emissionTmp[path_[i]]*(1.0 - missCopyProb), // probability of the same
                                   emissionTmp[(size_t)(1 -path_[i])] * missCopyProb }); // probability of differ

    normalizeBySum(sameDiffDist);

    if ( sampleIndexGivenProp( missCopyRg_, sameDiffDist) == 1 ) {

      hap_.push_back( 1 - path_[i] ); // differ
      siteOfOneMissCopyOne[i] += 1.0;

    } else {

      hap_.push_back( path_[i] ); // same

    }

  }

  assert ( hap_.size() == nLoci_ );

}


void UpdateSingleHap::sampleHapIndependently( std::vector <double> &plaf ){

  assert( hap_.size() == 0 );

  size_t plafIndex = segmentStartIndex_;

  for ( size_t i = 0; i < nLoci_; i++) {

    double tmpMax = max_value ( std::vector <double> ( {llk0_[i], llk1_[i]} ) ) ;
    std::vector <double> tmpDist ( {exp(llk0_[i] - tmpMax) * (1.0-plaf[plafIndex]),
                               exp(llk1_[i] - tmpMax) * plaf[plafIndex] } );
    (void)normalizeBySum(tmpDist);
    hap_.push_back ( (double)sampleIndexGivenProp(recombRg_, tmpDist) );
    plafIndex++;

  }

  assert ( hap_.size() == nLoci_ );

}


void UpdateSingleHap::updateLLK(){

  newLLK = std::vector <double> (nLoci_, 0.0);

  for ( size_t i = 0; i < nLoci_; i++){
    if ( hap_[i] == 0){

      newLLK[i] = llk0_[i];

    } else if (hap_[i] == 1){

      newLLK[i] = llk1_[i];

    } else {

      throw ShouldNotBeCalled();

    }

  }

}


UpdatePairHap::UpdatePairHap( std::vector <double> &refCount,
                              std::vector <double> &altCount,
                              std::vector <double> &plaf,
                              std::vector <double> &expectedWsaf,
                              std::vector <double> &proportion,
                              std::vector < std::vector <double> > &haplotypes,
                              RandomGenerator* rg,
                              size_t segmentStartIndex,
                              size_t nLoci,
                              std::shared_ptr<Panel> panel, double missCopyProb, double scalingFactor,
                              bool forbidCopyFromSame,
                              size_t strainIndex1,
                              size_t strainIndex2 ): UpdateHap(refCount,
                                                               altCount,
                                                               plaf,
                                                               expectedWsaf,
                                                               proportion,
                                                               haplotypes,
                                                               rg,
                                                               segmentStartIndex,
                                                               nLoci,
                                                               panel,
                                                               missCopyProb,
                                                               scalingFactor) {
    strainIndex1_ = strainIndex1;
    strainIndex2_ = strainIndex2;
    forbidCopyFromSame_ = forbidCopyFromSame;
    siteOfTwoSwitchOne = std::vector <double> (nLoci);
    siteOfTwoMissCopyOne = std::vector <double> (nLoci);
    siteOfTwoSwitchTwo = std::vector <double> (nLoci);
    siteOfTwoMissCopyTwo = std::vector <double> (nLoci);

}


void UpdatePairHap::core(std::vector <double> &refCount,
                           std::vector <double> &altCount,
                           std::vector <double> &plaf,
                           std::vector <double> &expectedWsaf,
                           std::vector <double> &proportion,
                           std::vector < std::vector <double> > &haplotypes){

  calcExpectedWsaf( expectedWsaf, proportion, haplotypes);
  calcHapLLKs(refCount, altCount);

  if ( panel_ != NULL ){

    buildEmission(missCopyProb_);
    calcFwdProbs(forbidCopyFromSame_);
    samplePaths();
    addMissCopying(missCopyProb_);

  } else {

    sampleHapIndependently( plaf );

  }

  updateLLK();

}


void UpdatePairHap:: calcExpectedWsaf( std::vector <double> & expectedWsaf, std::vector <double> &proportion, std::vector < std::vector <double> > &haplotypes) {
  expectedWsaf00_ = std::vector <double> (expectedWsaf.begin()+segmentStartIndex_, expectedWsaf.begin()+(segmentStartIndex_+nLoci_));

  size_t hapIndex = segmentStartIndex_;
  for ( size_t i = 0; i < expectedWsaf00_.size(); i++ ){

    expectedWsaf00_[i] -= (proportion[strainIndex1_] * haplotypes[hapIndex][strainIndex1_] + proportion[strainIndex2_] * haplotypes[hapIndex][strainIndex2_]);
    assert (expectedWsaf00_[i] >= 0 );
    assert (expectedWsaf00_[i] < 1 );
    hapIndex++;

  }

  expectedWsaf10_ = expectedWsaf00_;
  for ( size_t i = 0; i < expectedWsaf10_.size(); i++ ) {

    expectedWsaf10_[i] += proportion[strainIndex1_] ;

  }

  expectedWsaf01_ = expectedWsaf00_;
  for ( size_t i = 0; i < expectedWsaf01_.size(); i++ ) {

    expectedWsaf01_[i] += proportion[strainIndex2_] ;

  }

  expectedWsaf11_ = expectedWsaf00_;
  for ( size_t i = 0; i < expectedWsaf11_.size(); i++ ) {

    expectedWsaf11_[i] += (proportion[strainIndex1_] + proportion[strainIndex2_]);

  }

}


void UpdatePairHap:: calcHapLLKs( std::vector <double> &refCount, std::vector <double> &altCount) {

  llk00_ = calcLLKs( refCount, altCount, expectedWsaf00_, segmentStartIndex_, nLoci_, scalingFactor() );
  llk10_ = calcLLKs( refCount, altCount, expectedWsaf10_, segmentStartIndex_, nLoci_, scalingFactor() );
  llk01_ = calcLLKs( refCount, altCount, expectedWsaf01_, segmentStartIndex_, nLoci_, scalingFactor() );
  llk11_ = calcLLKs( refCount, altCount, expectedWsaf11_, segmentStartIndex_, nLoci_, scalingFactor() );
  assert( llk00_.size() == nLoci_ );
  assert( llk10_.size() == nLoci_ );
  assert( llk01_.size() == nLoci_ );
  assert( llk11_.size() == nLoci_ );

}


void UpdatePairHap:: buildEmission( double missCopyProb ) {

  std::vector <double> noMissProb (nLoci_, log(1.0 - missCopyProb));
  std::vector <double> missProb (nLoci_, log( missCopyProb ));
  std::vector <double> noNo = vecSum(noMissProb, noMissProb);
  std::vector <double> misMis = vecSum(missProb, missProb);
  std::vector <double> misNo = vecSum(noMissProb, missProb);

  std::vector <double> tmp_00_1 = vecSum(llk00_, noNo);
  std::vector <double> tmp_00_2 = vecSum(llk10_, misNo);
  std::vector <double> tmp_00_3 = vecSum(llk01_, misNo);
  std::vector <double> tmp_00_4 = vecSum(llk11_, misMis);

  std::vector <double> tmp_01_1 = vecSum(llk01_, noNo);
  std::vector <double> tmp_01_2 = vecSum(llk00_, misNo);
  std::vector <double> tmp_01_3 = vecSum(llk11_, misNo);
  std::vector <double> tmp_01_4 = vecSum(llk10_, misMis);

  std::vector <double> tmp_10_1 = vecSum(llk10_, noNo);
  std::vector <double> tmp_10_2 = vecSum(llk00_, misNo);
  std::vector <double> tmp_10_3 = vecSum(llk11_, misNo);
  std::vector <double> tmp_10_4 = vecSum(llk01_, misMis);

  std::vector <double> tmp_11_1 = vecSum(llk11_, noNo);
  std::vector <double> tmp_11_2 = vecSum(llk10_, misNo);
  std::vector <double> tmp_11_3 = vecSum(llk01_, misNo);
  std::vector <double> tmp_11_4 = vecSum(llk00_, misMis);

  assert(emission_.size() == 0 );

  for ( size_t i = 0; i < nLoci_; i++) {

    std::vector <double> tmp ({tmp_00_1[i], tmp_00_2[i], tmp_00_3[i], tmp_00_4[i],
                          tmp_01_1[i], tmp_01_2[i], tmp_01_3[i], tmp_01_4[i],
                          tmp_10_1[i], tmp_10_2[i], tmp_10_3[i], tmp_10_4[i],
                          tmp_11_1[i], tmp_11_2[i], tmp_11_3[i], tmp_11_4[i]});

    double tmaxTmp = max_value(tmp);

    std::vector <double> emissRow ({exp(tmp_00_1[i] - tmaxTmp) + exp(tmp_00_2[i] - tmaxTmp) + exp(tmp_00_3[i] - tmaxTmp) + exp(tmp_00_4[i] - tmaxTmp),
                               exp(tmp_01_1[i] - tmaxTmp) + exp(tmp_01_2[i] - tmaxTmp) + exp(tmp_01_3[i] - tmaxTmp) + exp(tmp_01_4[i] - tmaxTmp),
                               exp(tmp_10_1[i] - tmaxTmp) + exp(tmp_10_2[i] - tmaxTmp) + exp(tmp_10_3[i] - tmaxTmp) + exp(tmp_10_4[i] - tmaxTmp),
                               exp(tmp_11_1[i] - tmaxTmp) + exp(tmp_11_2[i] - tmaxTmp) + exp(tmp_11_3[i] - tmaxTmp) + exp(tmp_11_4[i] - tmaxTmp)});

    emission_.push_back(emissRow);

  }

    assert(emission_.size() == nLoci_ );

}


std::vector <double> UpdatePairHap::computeRowMarginalDist( std::vector < std::vector < double > > & probDist ){ // Sum of Rows

  std::vector <double> marginalDist (probDist.size(), 0.0);

  for ( size_t i = 0; i < probDist.size(); i++ ){

    marginalDist[i] = sumOfVec(probDist[i]);

  }

  return marginalDist;

}


std::vector <double> UpdatePairHap::computeColMarginalDist( std::vector < std::vector < double > > & probDist ){ // Sum of Cols

  std::vector <double> marginalDist (probDist.size(), 0.0);

  for ( size_t coli = 0; coli < probDist[0].size(); coli++ ) {

    for ( size_t rowi = 0; rowi < probDist.size(); rowi++ ) {

      marginalDist[coli] += probDist[rowi][coli];

    }

  }

  return marginalDist;

}


void UpdatePairHap:: calcFwdProbs( bool forbidCopyFromSame ){

  size_t hapIndex = segmentStartIndex_;
  assert ( fwdProbs_.size() == 0 );
  std::vector < std::vector < double > > fwd1st;

  for ( size_t i = 0 ; i < nPanel_; i++) { // Row of the matrix

    size_t rowObs = (size_t)panel_->content_[0][i];
    std::vector <double> fwd1stRow (nPanel_, 0.0);

    for ( size_t ii = 0 ; ii < nPanel_; ii++) { // Column of the matrix

      if ( forbidCopyFromSame && i == ii ) continue;

      size_t colObs = (size_t)panel_->content_[hapIndex][ii];
      size_t obs = rowObs*2 + colObs;
      fwd1stRow[ii] = emission_[0][obs];

    }

    fwd1st.push_back(fwd1stRow);

  }

  normalizeBySumMat(fwd1st);
  fwdProbs_.push_back(fwd1st);

  for ( size_t j = 1; j < nLoci_; j++ ) {

    double recRec = panel_->pRecRec_[hapIndex];
    double recNorec = panel_->pRecNoRec_[hapIndex];
    double norecNorec = panel_->pNoRecNoRec_[hapIndex];
    hapIndex++;

    std::vector <double> marginalOfRows = computeRowMarginalDist( fwdProbs_.back() );
    std::vector <double> marginalOfCols = computeColMarginalDist( fwdProbs_.back() );

    std::vector < std::vector < double > > fwdTmp;
    for ( size_t i = 0 ; i < nPanel_; i++) {

      size_t rowObs = (size_t)panel_->content_[hapIndex][i];
      std::vector <double> fwdTmpRow (nPanel_, 0.0);
      for ( size_t ii = 0 ; ii < nPanel_; ii++){
          if ( forbidCopyFromSame && i == ii ) continue;

          size_t colObs = (size_t)panel_->content_[hapIndex][ii];
          size_t obs = rowObs*2 + colObs;
          fwdTmpRow[ii] = emission_[j][obs] * (sumOfMat(fwdProbs_.back())*recRec +
                                                     fwdProbs_.back()[i][ii]*norecNorec+
                                                     recNorec * ( marginalOfRows[ii]+marginalOfCols[i] ) );

      }

      fwdTmp.push_back(fwdTmpRow);

    }

    normalizeBySumMat(fwdTmp);
    fwdProbs_.push_back(fwdTmp);

  }

}


std::vector <size_t> UpdatePairHap::sampleMatrixIndex( std::vector < std::vector < double > > &probDist ) {

  size_t tmp = sampleIndexGivenProp ( recombLevel2Rg_, reshapeMatToVec(probDist));
  div_t divresult;
  divresult = div((int)tmp, (int)nPanel_);

  return std::vector <size_t> ({(size_t)divresult.quot, (size_t)divresult.rem});

}


void UpdatePairHap::samplePaths() {

  assert ( path1_.size() == 0 );
  assert ( path2_.size() == 0 );

  std::vector <size_t> tmpPath = sampleMatrixIndex(fwdProbs_[nLoci_-1]);
  size_t rowI = tmpPath[0];
  size_t colJ = tmpPath[1];
  size_t contentIndex = segmentStartIndex_ + nLoci_ - 1;

  path1_.push_back(panel_->content_[contentIndex][rowI]);
  path2_.push_back(panel_->content_[contentIndex][colJ]);

  for ( size_t j = (nLoci_ - 1) ; j > 0; j-- ) {

    --contentIndex;
    double recRec = panel_->pRecRec_[contentIndex];
    double recNorec = panel_->pRecNoRec_[contentIndex];
    double norecNorec = panel_->pNoRecNoRec_[contentIndex];

    size_t previous_site = j - 1;
    std::vector < std::vector < double > > previousDist = fwdProbs_[previous_site];
    double previousProbij =previousDist[rowI][colJ];

    std::vector <double> rowIdist = previousDist[rowI];
    double tmpRowSum = sumOfVec(rowIdist);

    std::vector <double> colJdist;
    for ( auto const& array: previousDist ){

      colJdist.push_back( array[colJ] );

    }

    assert(nPanel_ == colJdist.size());

    double tmpColSum = sumOfVec(colJdist);

    std::vector <double> weightOfFourCases ({ recRec     * sumOfMat(previousDist),           // recombination happened on both strains
                                         recNorec   * tmpRowSum,  // first strain no recombine, second strain recombine
                                         recNorec   * tmpColSum,  // first strain recombine, second strain no recombine
                                         norecNorec * previousProbij }); // no recombine on either strain

    normalizeBySum(weightOfFourCases);

    size_t tmpCase = sampleIndexGivenProp( recombRg_, weightOfFourCases );

    if ( tmpCase == (size_t)0 ){ // switching both strains

      siteOfTwoSwitchTwo[j] += 1.0;
      tmpPath = sampleMatrixIndex(previousDist);
      rowI = tmpPath[0];
      colJ = tmpPath[1];
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else if ( tmpCase == 1 ) { // switching second strain

      siteOfTwoSwitchOne[j] += 0.5;
      normalizeBySum(rowIdist);
      colJ = sampleIndexGivenProp( recombLevel2Rg_, rowIdist );
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else if ( tmpCase == (size_t)2 ){ // switching first strain

      siteOfTwoSwitchOne[j] += 0.5;
      normalizeBySum(colJdist);
      rowI = sampleIndexGivenProp( recombLevel2Rg_, colJdist );
      colJ = colJ;
      //assert (rowI != colJ); // OFF, as by default, allow copying the same strain

    } else if ( tmpCase == 3 ) { // no switching

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


void UpdatePairHap::addMissCopying( double missCopyProb ) {

  assert( hap1_.size() == 0 );
  assert( hap2_.size() == 0 );

  for ( size_t i = 0; i < nLoci_; i++){

    double tmpMax = max_value ( std::vector <double>({llk00_[i], llk01_[i], llk10_[i], llk11_[i]}));

    std::vector <double> emissionTmp ({exp(llk00_[i]-tmpMax), exp(llk01_[i]-tmpMax), exp(llk10_[i]-tmpMax), exp(llk11_[i]-tmpMax)});

    std::vector <double> casesDist ( { emissionTmp[(size_t)(2*path1_[i]     +path2_[i])]     * (1.0 - missCopyProb) * (1.0 - missCopyProb), // probability of both same
                                  emissionTmp[(size_t)(2*path1_[i]     +(1-path2_[i]))] * (1.0 - missCopyProb) * missCopyProb,         // probability of same1diff2
                                  emissionTmp[(size_t)(2*(1 -path1_[i])+path2_[i])]     * missCopyProb * (1.0 - missCopyProb),         // probability of same2diff1
                                  emissionTmp[(size_t)(2*(1 -path1_[i])+(1-path2_[i]))] * missCopyProb * missCopyProb });              // probability of both differ
    normalizeBySum(casesDist);

    size_t tmpCase = sampleIndexGivenProp( missCopyRg_, casesDist );

    if ( tmpCase == 0 ){

      hap1_.push_back( path1_[i] );
      hap2_.push_back( path2_[i] );

    } else if ( tmpCase == 1 ){

      siteOfTwoMissCopyOne[i] += 0.5;
      hap1_.push_back( path1_[i] );
      hap2_.push_back( 1.0 - path2_[i] );

    } else if ( tmpCase == 2 ){

      siteOfTwoMissCopyOne[i] += 0.5;
      hap1_.push_back( 1.0 - path1_[i] );
      hap2_.push_back( path2_[i] );

    } else if ( tmpCase == 3 ){

      siteOfTwoMissCopyTwo[i] += 1.0;
      hap1_.push_back( 1.0 - path1_[i] );
      hap2_.push_back( 1.0 - path2_[i] );

    } else {

      throw ShouldNotBeCalled();

    }

  }

  assert ( hap1_.size() == nLoci_ );
  assert ( hap2_.size() == nLoci_ );

}


void UpdatePairHap::sampleHapIndependently(std::vector <double> &plaf){

  assert( hap1_.size() == 0 );
  assert( hap2_.size() == 0 );

  size_t plafIndex = segmentStartIndex_;

  for ( size_t i = 0; i < nLoci_; i++){

    double tmpMax = max_value ( std::vector <double> ( {llk00_[i], llk01_[i], llk10_[i], llk11_[i]} ) );

    std::vector <double> tmpDist ( {exp(llk00_[i] - tmpMax) * (1.0-plaf[plafIndex]) * (1.0-plaf[plafIndex]),
                               exp(llk01_[i] - tmpMax) * (1.0-plaf[plafIndex]) * plaf[plafIndex],
                               exp(llk10_[i] - tmpMax) * (1.0-plaf[plafIndex]) * plaf[plafIndex],
                               exp(llk11_[i] - tmpMax) * plaf[plafIndex] * plaf[plafIndex] } );

    normalizeBySum(tmpDist);

    size_t tmpCase = sampleIndexGivenProp( recombRg_, tmpDist );

    if ( tmpCase == 0 ){

      hap1_.push_back( 0.0 );
      hap2_.push_back( 0.0 );

    } else if ( tmpCase == 1 ){

      hap1_.push_back( 0.0 );
      hap2_.push_back( 1.0 );

    } else if ( tmpCase == 2 ){

      hap1_.push_back( 1.0 );
      hap2_.push_back( 0.0 );

    } else if ( tmpCase == 3 ){

      hap1_.push_back( 1.0 );
      hap2_.push_back( 1.0 );

    } else {

      throw ShouldNotBeCalled();

    }

    plafIndex++;

  }

  assert ( hap1_.size() == nLoci_ );
  assert ( hap2_.size() == nLoci_ );

}


void UpdatePairHap::updateLLK(){

  newLLK = std::vector <double> (nLoci_, 0.0);

  for ( size_t i = 0; i < nLoci_; i++) {

    if ( hap1_[i] == 0 && hap2_[i] == 0 ) {

      newLLK[i] = llk00_[i];

    } else if (hap1_[i] == 0 && hap2_[i] == 1){

      newLLK[i] = llk01_[i];

    } else if (hap1_[i] == 1 && hap2_[i] == 0){

      newLLK[i] = llk10_[i];

    } else if (hap1_[i] == 1 && hap2_[i] == 1){

      newLLK[i] = llk11_[i];

    } else {

      throw ShouldNotBeCalled();

    }

  }

}
