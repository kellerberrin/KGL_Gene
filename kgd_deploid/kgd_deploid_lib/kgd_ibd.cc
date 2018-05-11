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

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "kgd_ibd.h"


namespace kgd = kellerberrin::deploid;



void kgd::IBDconfiguration::buildIBDconfiguration(size_t k) {

  this->setKstrain(k);
  this->enumerateOp();
  this->makePairList();
  this->makePairToEmission();
  this->findUniqueState();
  this->findEffectiveK();

}



void kgd::IBDconfiguration::enumerateOp() {
  //#For each configuration, identify which pairs are IBD

  this->op = enumerateBinaryMatrixOfK(nchoose2(this->kStrain()));

}


void kgd::IBDconfiguration::makePairList() {
  //#Make a map of pairs to pair value

  assert(pairList.size() == 0);

  for (size_t i = 0; i < this->kStrain(); i++) { // 0-indexed

    for (size_t j = i + 1; j < this->kStrain(); j++) {

      pairList.push_back(std::vector<size_t>({i, j}));

    }

  }

  assert(pairList.size() == (size_t) nchoose2(this->kStrain()));

}


void kgd::IBDconfiguration::makePairToEmission() {

  assert(pairToEmission.size() == 0);
  //prs2ems<-array(0, c(nrow(op), k));
  //for ( size_t i = 0; i < this->op.size(); i++ ){

  for (std::vector<int> tmpOp : op) {

    std::vector<int> tmpRow = makeTmpRow();
    //for ( size_t i = 0 ; i <(*opIt).size(); i++){
    //cout << (*opIt)[i] << " ";
    //}
    //cout<<endl;

    std::vector<size_t> ii = findWhichIsOne(tmpOp);
    //ii <- which(op[rowI,]==1);
    //cout << ii.size()<<endl;
    //cout << "##############" <<endl;

    if (ii.size() > 0) {

      std::vector<std::vector<size_t> > tmpIBDPairs;

      for (size_t j = 0; j < ii.size(); j++) {
        //cout << "j = " << j <<" ii[j] = "<<ii[j]<<endl;
        tmpIBDPairs.push_back(pairList[ii[j]]);

        //cout << tmpIBDPairs.back()[0]<< " "<<tmpIBDPairs.back()[1]<<endl;
      }

      int tmpIndex = (tmpIBDPairs.size() - 1);

      while (tmpIndex >= 0) {

        //cout << "replacing element "<< tmpIBDPairs[tmpIndex][0] << " by element " << tmpIBDPairs[tmpIndex][1] <<endl;
        tmpRow[tmpIBDPairs[tmpIndex][0]] = tmpRow[tmpIBDPairs[tmpIndex][1]];
        tmpIndex--;

      }

    }

    pairToEmission.push_back(tmpRow);
    //for (size_t i = 0 ; i < pairToEmission.back().size(); i++){
    //cout << pairToEmission.back()[i]<<" ";
    //}
    //cout <<endl;
  }
  //for (size_t i = 0; i < pairToEmission.size(); i++){
  //for (size_t ii = 0 ; ii < pairToEmission.back().size(); ii++){
  //cout << pairToEmission[i][ii]<<" ";
  //}
  //cout <<endl;
  //}

}


std::vector<int> kgd::IBDconfiguration::makeTmpRow() {

  std::vector<int> ret(this->kStrain());

  for (size_t i = 0; i < ret.size(); i++) {
    ret[i] = (int) i;
  }

  return ret;

}


std::vector<size_t> kgd::IBDconfiguration::findWhichIsOne(std::vector<int> tmpOp) {

  std::vector<size_t> ret;

  for (size_t i = 0; i < tmpOp.size(); i++) {

    if (tmpOp[i] == 1) {

      ret.push_back(i);

    }

  }

  return ret;

}


void kgd::IBDconfiguration::findUniqueState() {

  assert (states.size() == 0);
  //states.push_back(pairToEmission[0]);
  //for (size_t i = 1; i < this->pairToEmission.size(); i++){
  //bool aNewState = true;
  //for ( vector<int> state : states){
  //if ( twoVectorsAreSame(state, this->pairToEmission[i]) ){
  //aNewState = false;
  //break;
  //}
  //}
  //if ( aNewState ){
  //states.push_back(this->pairToEmission[i]);
  //}
  //}

  states = unique(this->pairToEmission);

  //for ( vector<int> state : states){
  //for (int i : state){
  //cout << i <<" ";
  //}
  //cout <<endl;
  //}

}


void kgd::IBDconfiguration::findEffectiveK() {

  assert(effectiveK.size() == 0);

  for (std::vector<int> state : states) {

    std::set<int> tmpSet(state.begin(), state.end());
    //cout << tmpSet.size() <<endl;
    effectiveK.push_back(tmpSet.size());

  }

  assert(effectiveK.size() == states.size());

}


std::vector<std::string> kgd::IBDconfiguration::getIBDconfigureHeader() {

  std::vector<std::string> ret;

  for (size_t i = 0; i < this->states.size(); i++) {

    std::string tmp;

    for (size_t j = 0; j < this->states[i].size(); j++) {

      std::stringstream tmp_ss;
      tmp_ss << this->states[i][j];
      tmp += tmp_ss.str() + ((j < (this->states[i].size() - 1)) ? "-" : "");

    }

    ret.push_back(tmp);

  }

  return ret;

}



void kgd::Hprior::buildHprior(size_t kStrain, std::vector<double> &plaf) {

  ibdConfig.buildIBDconfiguration(kStrain);

  this->effectiveK = ibdConfig.effectiveK;
  this->nState_ = 0;
  this->plaf_ = plaf;
  this->setKstrain(kStrain);
  this->setnLoci(this->plaf_.size());

  std::vector<std::vector<int> > hSetBase = enumerateBinaryMatrixOfK(this->kStrain());
  size_t stateI = 0;

  for (std::vector<int> state : ibdConfig.states) {

    std::set<int> stateUnique(state.begin(), state.end());

    assert(stateUnique.size() == effectiveK[stateI]);

    std::vector<std::vector<int> > hSetBaseTmp = hSetBase;

    for (size_t j = 0; j < (size_t) this->kStrain(); j++) {

      for (size_t i = 0; i < hSetBase.size(); i++) {

        hSetBaseTmp[i][j] = hSetBase[i][state[j]];

      }

    }

    std::vector<std::vector<int> > hSetBaseTmpUnique = unique(hSetBaseTmp); // uu

    size_t sizeOfhSetBaseTmpUnique = hSetBaseTmpUnique.size();

    stateIdxFreq.push_back(sizeOfhSetBaseTmpUnique);
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
        //cout << (1.0-plaf_[site]) << " " << tmpDiff << " " <<pow((1.0-plaf_[site]),(double)tmpDiff) << endl;

        hPriorTmp[site] = pow(plaf_[site], (double) tmpSum) * pow((1.0 - plaf_[site]), (double) tmpDiff);

      }

      priorProb.push_back(hPriorTmp);
      hSet.push_back(hSetBaseTmpUnique[i]);

      nState_++;
      stateIdx.push_back(stateI);

    }

    stateI++;

  }

}


void kgd::Hprior::transposePriorProbs() {

  assert(priorProbTrans.size() == 0);

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> priorProbTransTmp(nState());

    for (size_t j = 0; j < nState(); j++) {

      priorProbTransTmp[j] = priorProb[j][i];
      //cout << priorProb[j][i] << " ";
    }

    priorProbTrans.push_back(priorProbTransTmp);
    //cout << endl;
  }

}


std::vector<std::string> kgd::Hprior::getIBDconfigureHeader() {

  return this->ibdConfig.getIBDconfigureHeader();

}


std::vector<std::vector<int> > kgd::unique(std::vector<std::vector<int> > &mat) {

  std::vector<std::vector<int> > ret;

  ret.push_back(mat[0]);

  for (size_t i = 1; i < mat.size(); i++) {

    bool aNewState = true;

    for (std::vector<int> state : ret) {

      if (twoVectorsAreSame(state, mat[i])) {
        aNewState = false;
        break;
      }

    }

    if (aNewState) {
      ret.push_back(mat[i]);
    }

  }

  return ret;

}


std::vector<std::vector<int> > kgd::enumerateBinaryMatrixOfK(size_t k) {
  // This function enumerate all possible binary combinations of k elements

  int ksq = pow(2, k);

  std::vector<std::vector<int> > ret;

  for (int i = 0; i < ksq; i++) {

    ret.push_back(convertIntToBinary(i, k));

  }

  return ret;

}


std::vector<int> kgd::convertIntToBinary(int x, size_t len) {

  std::vector<int> ret(len);

  size_t idx = 0;

  while (x) {

    ret[idx] = (x & 1) ? 1 : 0;
    idx++;
    //cout << "x " <<x<< " idx "<<idx<<" len "<< len<<endl;

    if (idx > len) {

      throw OutOfVectorSize();
    }

    x >>= 1;

  }

  std::reverse(ret.begin(), ret.end());
  //for (size_t i = 0; i < ret.size(); i++){
  //cout << ret[i] << " ";
  //}
  //cout<<endl;
  return ret;
}


int kgd::nchoose2(int n) {

  if (n < 2) {

    throw InvalidInput("Input must be at least 2!");

  }

  int ret = n * (n - 1) / 2;
  return ret;

}


bool kgd::twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2) {

  if (vec1.size() != vec2.size()) {

    throw InvalidInput("Input vectors have different length!");

  }

  bool ret = true;

  for (size_t i = 0; i < vec1.size(); i++) {

    if (vec1[i] != vec2[i]) {

      ret = false;
      break;

    }

  }

  return ret;

}



void kgd::IBDpath::init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator) {

  this->ibdRg_ = randomGenerator;
  this->setNLoci(dEploidIO.nLoci());
  this->setKstrain(dEploidIO.kStrain());
  this->setTheta(1.0 / (double) kStrain());

  this->IBDpathChangeAt = std::vector<double>(this->nLoci());

  // compute likelihood surface
  this->makeLlkSurf(dEploidIO.altCount_, dEploidIO.refCount_);

  // initialize haplotype prior
  this->hprior.buildHprior(kStrain(), dEploidIO.plaf_);
  this->hprior.transposePriorProbs();

  this->makeIbdTransProbs();

  // initialize fm
  this->fSumState = std::vector<double>(this->hprior.nPattern());

  // initialize ibdConfigurePath
  this->ibdConfigurePath = std::vector<size_t>(this->nLoci());

  // initialize recombination probabilities;
  this->ibdRecombProbs = IBDrecombProbs(dEploidIO.position_, dEploidIO.nLoci());

  this->ibdRecombProbs.computeRecombProbs(dEploidIO.averageCentimorganDistance(),
                                          dEploidIO.parameterG(),
                                          dEploidIO.useConstRecomb(),
                                          dEploidIO.constRecombProb());

  this->currentIBDpathChangeAt = std::vector<double>(this->nLoci());

  this->computeUniqueEffectiveKCount();

}


void kgd::IBDpath::ibdSamplePath(std::vector<double> statePrior) {

  int lociIdx = this->nLoci() - 1;

  std::vector<double> tmpProp = fm[lociIdx];
  (void) normalizeBySum(tmpProp);

  ibdConfigurePath[lociIdx] = sampleIndexGivenProp(this->ibdRg_, tmpProp);

  assert(this->fm.size() == nLoci());

  while (lociIdx > 0) {

    lociIdx--;

    std::vector<double> vNoRecomb = vecProd(this->ibdTransProbs[this->hprior.stateIdx[ibdConfigurePath[lociIdx + 1]]],
                                       fm[lociIdx]);
    assert(vNoRecomb.size() == this->hprior.nState());

    std::vector<double> vRecomb = fm[lociIdx];

    assert(vRecomb.size() == this->hprior.nState());

    std::vector<double> prop(this->hprior.nState());

    for (size_t i = 0; i < prop.size(); i++) {

      prop[i] = vNoRecomb[i] * this->ibdRecombProbs.pNoRec_[lociIdx] +
                vRecomb[i] * this->ibdRecombProbs.pRec_[lociIdx] * statePrior[ibdConfigurePath[lociIdx + 1]];

    }

    tmpProp = prop;
    normalizeBySum(tmpProp);
    ibdConfigurePath[lociIdx] = sampleIndexGivenProp(this->ibdRg_, tmpProp);

    assert(ibdConfigurePath[lociIdx] < this->hprior.nState());

  }

}


std::vector<size_t> kgd::IBDpath::findWhichIsSomething(std::vector<size_t> tmpOp, size_t something) {

  std::vector<size_t> ret;

  for (size_t i = 0; i < tmpOp.size(); i++) {

    if (tmpOp[i] == something) {
      ret.push_back(i);
    }

  }

  return ret;

}


void kgd::IBDpath::buildPathProbabilityForPainting(std::vector<double> proportion) {

  //vector <double> effectiveKPrior = this->computeEffectiveKPrior(this->theta());
  std::vector<double> effectiveKPrior = std::vector<double>(this->hprior.nPattern(), 1.0 / this->hprior.nPattern());
  std::vector<double> statePrior = this->computeStatePrior(effectiveKPrior);

  // First building the path likelihood
  this->computeIbdPathFwdProb(proportion, statePrior);

  // Reshape Fwd
  std::vector<std::vector<double>> reshapedFwd = reshapeProbs(this->fm);

  this->computeIbdPathBwdProb(proportion, effectiveKPrior, statePrior);
  // Reshape Bwd
  std::vector<std::vector<double>> reshapedBwd = reshapeProbs(this->bwd);

  // Combine Fwd Bwd
  this->combineFwdBwd(reshapedFwd, reshapedBwd);

}


void kgd::IBDpath::computeIbdPathBwdProb(std::vector<double> proportion, std::vector<double> effectiveKPrior, std::vector<double> statePrior) {
  //# assuming each ibd state has equal probabilities, transform it into ibd configurations
  //dout << " start building ibd bwd "<< endl;
  std::vector<double> tmp = std::vector<double>(hprior.stateIdxFreq.size());

  assert(effectiveKPrior.size() == hprior.stateIdxFreq.size());

  for (size_t i = 0; i < tmp.size(); i++) {

    tmp[i] = effectiveKPrior[i] / (double) hprior.stateIdxFreq[i];

  }

  std::vector<double> tmpBw = std::vector<double>(hprior.nState());

  for (size_t j = 0; j < tmpBw.size(); j++) {

    for (size_t i = 0; i < tmp.size(); i++) {

      tmpBw[j] += tmp[i] * ibdTransProbs[i][j];

    }

  }

  this->bwd.push_back(tmpBw);

  for (size_t rev_siteI = 1; rev_siteI < this->nLoci(); rev_siteI++) {

    size_t siteI = this->nLoci() - rev_siteI;

    std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, siteI);

    //vector<double> lk = vector <double> (hprior.nState(), 1.0);
    std::vector<double> bSumState = std::vector<double>(hprior.nPattern());
    for (size_t i = 0; i < bSumState.size(); i++) {

      for (size_t j = 0; j < hprior.nState(); j++) {

        bSumState[i] += ibdTransProbs[i][j] * this->bwd.back()[j];

      }

    }

    std::vector<double> vNoRecomb(hprior.nState());
    for (size_t i = 0; i < hprior.stateIdx.size(); i++) {

      vNoRecomb[i] = bSumState[hprior.stateIdx[i]];

    }

    for (size_t i = 0; i < hprior.nState(); i++) {

      tmpBw[i] = 0;

      for (size_t j = 0; j < lk.size(); j++) {

        tmpBw[i] += (lk[j] * bwd.back()[j]) * this->ibdRecombProbs.pRec_[siteI - 1];

      }

      tmpBw[i] *= statePrior[i];
      tmpBw[i] += lk[i] * (this->ibdRecombProbs.pNoRec_[siteI - 1]) * vNoRecomb[i];
      tmpBw[i] *= hprior.priorProb[i][siteI];

    }

    normalizeBySum(tmpBw);
    this->bwd.push_back(tmpBw);

  }

  reverse(bwd.begin(), bwd.end());

}


void kgd::IBDpath::computeIbdPathFwdProb(std::vector<double> proportion, std::vector<double> statePrior) {

  this->fm.clear();

  std::vector<double> vPrior = vecProd(statePrior, this->hprior.priorProbTrans[0]);

  std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, 0);

  this->updateFmAtSiteI(vPrior, lk);

  for (size_t siteI = 1; siteI < this->nLoci(); siteI++) {

    std::vector<double> vNoRec;

    for (size_t stateIdxTmp : hprior.stateIdx) {

      vNoRec.push_back(this->fSumState[stateIdxTmp]);

    }

    for (size_t i = 0; i < hprior.nState(); i++) {

      vPrior[i] = (vNoRec[i] * this->ibdRecombProbs.pNoRec_[siteI] + fSum * this->ibdRecombProbs.pRec_[siteI] * statePrior[i]) * hprior.priorProbTrans[siteI][i];

    }

    lk = computeLlkOfStatesAtSiteI(proportion, siteI);
    //cout << "lk = " ; for (double l :lk){cout << l << " ";}cout<<endl;

    this->updateFmAtSiteI(vPrior, lk);
    //for (double p : this->fm.back()){printf("%8.4f ", p);}cout<<endl;

  }

}


void kgd::IBDpath::updateFmAtSiteI(std::vector<double> &prior, std::vector<double> &llk) {

  std::vector<double> postAtSiteI = vecProd(prior, llk);
  //normalizeByMax(postAtSiteI);

  normalizeBySum(postAtSiteI);

  this->fm.push_back(postAtSiteI);
  this->fSum = sumOfVec(postAtSiteI);

  for (size_t i = 0; i < fSumState.size(); i++) {

    this->fSumState[i] = 0;

    for (size_t j = 0; j < hprior.nState(); j++) {

      this->fSumState[i] += ibdTransProbs[i][j] * postAtSiteI[j];

    }

  }

}


double kgd::IBDpath::bestPath(std::vector<double> proportion, double err) {

  double sumLLK = 0.0;

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> tmp;

    for (size_t j = 0; j < fm[i].size(); j++) {

      tmp.push_back(exp(log(fm[i][j]) + log(bwd[i][j])));

    }

    normalizeBySum(tmp);

    size_t indx = std::distance(tmp.begin(), std::max_element(tmp.begin(), tmp.end()));

    std::vector<int> hSetI = this->hprior.hSet[indx];

    double qs = 0;

    for (size_t j = 0; j < this->kStrain(); j++) {

      qs += (double) hSetI[j] * proportion[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;

    if ((qs > 0) & (qs < 1)) {

      sumLLK += logBetaPdf(qs2, this->llkSurf[i][0], this->llkSurf[i][1]);

    }

  }

  return sumLLK;

}


void kgd::IBDpath::combineFwdBwd(std::vector<std::vector<double>> &reshapedFwd, std::vector<std::vector<double>> &reshapedBwd) {

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> tmp;

    //cout << " site " << i << endl;

    for (size_t j = 0; j < reshapedFwd[i].size(); j++) {


      tmp.push_back(exp(log(reshapedFwd[i][j]) + log(reshapedBwd[i][j])));
      //tmp.push_back(exp(log(bwd[i][j])));
      //tmp.push_back(exp(log(fm[i][j])));
      //cout << "fwd = "<<fm[i][j]<<" "<< ", bwd = "<<bwd[i][j]<< ", fwdbwd = "<< tmp.back()<<endl;

    }

    normalizeBySum(tmp);
    fwdbwd.push_back(tmp);

  }
//post = exp(log(fmNomalized) + log(bwdNormalized))
}


void kgd::IBDpath::makeIbdTransProbs() {

  assert(this->ibdTransProbs.size() == 0);

  for (size_t i = 0; i < hprior.nPattern(); i++) {

    std::vector<double> transProbRow(hprior.nState());
    std::vector<size_t> wi = findWhichIsSomething(hprior.stateIdx, i);

    for (size_t wii : wi) {

      transProbRow[wii] = 1;

    }

    ibdTransProbs.push_back(transProbRow);

  }

}


std::vector<std::string> kgd::IBDpath::getIBDprobsHeader() {

  return this->hprior.getIBDconfigureHeader();

}


std::vector<std::vector<double> > kgd::IBDpath::reshapeProbs(std::vector<std::vector<double> > &probs) {

  assert(this->nLoci() == probs.size());

  std::vector<std::vector<double> > ret;

  for (size_t siteIndex = 0; siteIndex < this->nLoci(); siteIndex++) {

    size_t previousStateIdx = 0;
    std::vector<double> tmpRow;
    double cumProb = 0;

    for (size_t prob_ij = 0; prob_ij < probs[siteIndex].size(); prob_ij++) {

      cumProb += probs[siteIndex][prob_ij];

      if (previousStateIdx != this->hprior.stateIdx[prob_ij]) {
        cumProb -= probs[siteIndex][prob_ij];
        previousStateIdx++;
        tmpRow.push_back(cumProb);
        cumProb = probs[siteIndex][prob_ij];
      }
    }

    tmpRow.push_back(cumProb);
    normalizeBySum(tmpRow);
    ret.push_back(tmpRow);

  }

  return ret;

}


std::vector<double> kgd::IBDpath::computeEffectiveKPrior(double theta) {
  //#Calculate state prior given theta (theta is prob IBD)

  std::vector<double> pr0(this->kStrain());

  for (int i = 0; i < (int) pr0.size(); i++) {

    pr0[i] = binomialPdf(i, (int) (this->kStrain() - 1), theta);

  }

  std::vector<double> effectiveKPrior;

  for (size_t effectiveKtmp : this->hprior.effectiveK) {

    int effectiveKidx = effectiveKtmp - 1;
    assert(effectiveKidx >= 0);
    assert(effectiveKidx < (int) this->kStrain());
    effectiveKPrior.push_back(pr0[effectiveKidx] / uniqueEffectiveKCount[effectiveKidx]);

  }

  return effectiveKPrior;

}


std::vector<double> kgd::IBDpath::computeStatePrior(std::vector<double> effectiveKPrior) {

  std::vector<double> ret;

  for (size_t stateIdxTmp : this->hprior.stateIdx) {

    ret.push_back(effectiveKPrior[stateIdxTmp]);

  }

  return ret;

}


void kgd::IBDpath::computeAndUpdateTheta() {

  std::vector<size_t> obsState;

  size_t previousState = 0;
  size_t atSiteI = 0;

  for (size_t a : this->ibdConfigurePath) {

    if (a != previousState) {

      obsState.push_back(a);

    }

    if (this->hprior.stateIdx[a] != this->hprior.stateIdx[previousState]) {

      this->IBDpathChangeAt[atSiteI] += 1.0;
      this->currentIBDpathChangeAt[atSiteI] = 1.0;

    } else {

      this->currentIBDpathChangeAt[atSiteI] = 0.0;

    }

    previousState = a;
    atSiteI++;

  }

  size_t sumOfKeffStates = 0;
  size_t sccs = 0;

  for (size_t obs : obsState) {

    sumOfKeffStates += this->hprior.effectiveK[obs] - 1;
    sccs += this->kStrain() - this->hprior.effectiveK[obs];

  }
  //this->setTheta(rBeta(sccs+1.0, sumOfKeffStates+1.0, this->propRg_));

  this->setTheta(rBeta(sccs + 1.0, sumOfKeffStates + 1.0, this->ibdRg_));

}


void kgd::IBDpath::computeUniqueEffectiveKCount() {

  this->uniqueEffectiveKCount = std::vector<int>(this->kStrain());

  for (size_t effectiveKtmp : this->hprior.effectiveK) {

    int effectiveKidx = effectiveKtmp - 1;

    assert(effectiveKidx >= 0);

    this->uniqueEffectiveKCount[effectiveKidx]++;

  }

}


void kgd::IBDpath::makeLlkSurf(std::vector<double> altCount,
                               std::vector<double> refCount,
                               double scalingConst,
                               double err,
                               size_t gridSize) {

  double pGridSpacing = 1.0 / (double) (gridSize + 1);

  std::vector<double> pGrid;

  pGrid.push_back(pGridSpacing);

  for (size_t i = 1; i < gridSize; i++) {

    pGrid.push_back(pGrid.back() + pGridSpacing);

  }

  assert(pGrid.size() == gridSize);

  assert(llkSurf.size() == 0);

  for (size_t i = 0; i < altCount.size(); i++) {

    double alt = altCount[i];
    double ref = refCount[i];

    std::vector<double> ll;

    for (double unadjustedP : pGrid) {

      ll.push_back(calcLLK(ref, alt, unadjustedP, err, scalingConst));

    }

    double llmax = max_value(ll);
    std::vector<double> ln;

    for (double lltmp : ll) {

      ln.push_back(exp(lltmp - llmax));

    }

    double lnSum = sumOfVec(ln);
    for (size_t i = 0; i < ln.size(); i++) {

      ln[i] = ln[i] / lnSum;

    }

    std::vector<double> tmpVec1 = vecProd(ln, pGrid);
    double mn = sumOfVec(tmpVec1);
    std::vector<double> pGridSq = vecProd(pGrid, pGrid);
    std::vector<double> tmpVec2 = vecProd(ln, pGridSq);
    double vr = sumOfVec(tmpVec2) - mn * mn;

    double comm = (mn * (1.0 - mn) / vr - 1.0);
    llkSurf.push_back(std::vector<double>{mn * comm, (1 - mn) * comm});

  }

  assert(llkSurf.size() == this->nLoci());

}


std::vector<double> kgd::IBDpath::computeLlkOfStatesAtSiteI(std::vector<double> proportion, size_t siteI, double err) {

  std::vector<double> llks;

  for (std::vector<int> hSetI : this->hprior.hSet) {

    double qs = 0;

    for (size_t j = 0; j < this->kStrain(); j++) {

      qs += (double) hSetI[j] * proportion[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;
    //cout << qs2 << endl;

    llks.push_back(logBetaPdf(qs2, this->llkSurf[siteI][0], this->llkSurf[siteI][1]));

  }

  double maxllk = max_value(llks);
  std::vector<double> ret;

  for (double llk : llks) {

    double normalized = exp(llk - maxllk);
    //cout  << "llk-maxllk = "<<llk-maxllk<< " normalized " <<normalized<< endl;
    if (normalized == 0) {

      normalized = std::numeric_limits<double>::min();
      //cout << normalized<<endl;
    }

    ret.push_back(normalized);

  }

  return ret;
}

