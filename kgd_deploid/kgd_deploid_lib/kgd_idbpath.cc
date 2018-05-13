//
// Created by kellerberrin on 13/05/18.
//


#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <kgl_exec_env.h>
#include "kgd_ibdpath.h"



namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;




void kgd::IBDpath::init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator) {

  ibdRg_ = randomGenerator;
  setNLoci(dEploidIO.nLoci());
  setKstrain(dEploidIO.kStrain());
  setTheta(1.0 / (double) kStrain());

  IBDpathChangeAt = std::vector<double>(nLoci());

  // compute likelihood surface
  makeLlkSurf(dEploidIO.altCount_, dEploidIO.refCount_);

  // initialize haplotype prior
  hprior.buildHprior(kStrain(), dEploidIO.plaf_);
  hprior.transposePriorProbs();

  makeIbdTransProbs();

  // initialize fm
  fSumState = std::vector<double>(hprior.nPattern());

  // initialize ibdConfigurePath
  ibdConfigurePath = std::vector<size_t>(nLoci());

  // initialize recombination probabilities;
  ibdRecombProbs = IBDrecombProbs(dEploidIO.position_, dEploidIO.nLoci());

  ibdRecombProbs.computeRecombProbs(dEploidIO.averageCentimorganDistance(),
                                    dEploidIO.parameterG(),
                                    dEploidIO.useConstRecomb(),
                                    dEploidIO.constRecombProb());

  currentIBDpathChangeAt = std::vector<double>(nLoci());

  computeUniqueEffectiveKCount();

}


void kgd::IBDpath::ibdSamplePath(std::vector<double> statePrior) {

  int lociIdx = nLoci() - 1;

  std::vector<double> tmpProp = fm[lociIdx];
  (void) normalizeBySum(tmpProp);

  ibdConfigurePath[lociIdx] = sampleIndexGivenProp(ibdRg_, tmpProp);

  assert(fm.size() == nLoci());

  while (lociIdx > 0) {

    lociIdx--;

    std::vector<double> vNoRecomb = vecProd(ibdTransProbs[hprior.stateIdx[ibdConfigurePath[lociIdx + 1]]],
                                            fm[lociIdx]);
    assert(vNoRecomb.size() == hprior.nState());

    std::vector<double> vRecomb = fm[lociIdx];

    assert(vRecomb.size() == hprior.nState());

    std::vector<double> prop(hprior.nState());

    for (size_t i = 0; i < prop.size(); i++) {

      prop[i] = vNoRecomb[i] * ibdRecombProbs.pNoRec_[lociIdx] +
                vRecomb[i] * ibdRecombProbs.pRec_[lociIdx] * statePrior[ibdConfigurePath[lociIdx + 1]];

    }

    tmpProp = prop;
    normalizeBySum(tmpProp);
    ibdConfigurePath[lociIdx] = sampleIndexGivenProp(ibdRg_, tmpProp);

    assert(ibdConfigurePath[lociIdx] < hprior.nState());

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

  //vector <double> effectiveKPrior = computeEffectiveKPrior(theta());
  std::vector<double> effectiveKPrior = std::vector<double>(hprior.nPattern(), 1.0 / hprior.nPattern());
  std::vector<double> statePrior = computeStatePrior(effectiveKPrior);

  // First building the path likelihood
  computeIbdPathFwdProb(proportion, statePrior);

  // Reshape Fwd
  std::vector<std::vector<double>> reshapedFwd = reshapeProbs(fm);

  computeIbdPathBwdProb(proportion, effectiveKPrior, statePrior);
  // Reshape Bwd
  std::vector<std::vector<double>> reshapedBwd = reshapeProbs(bwd);

  // Combine Fwd Bwd
  combineFwdBwd(reshapedFwd, reshapedBwd);

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

  bwd.push_back(tmpBw);

  for (size_t rev_siteI = 1; rev_siteI < nLoci(); rev_siteI++) {

    size_t siteI = nLoci() - rev_siteI;

    std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, siteI);

    //vector<double> lk = vector <double> (hprior.nState(), 1.0);
    std::vector<double> bSumState = std::vector<double>(hprior.nPattern());
    for (size_t i = 0; i < bSumState.size(); i++) {

      for (size_t j = 0; j < hprior.nState(); j++) {

        bSumState[i] += ibdTransProbs[i][j] * bwd.back()[j];

      }

    }

    std::vector<double> vNoRecomb(hprior.nState());
    for (size_t i = 0; i < hprior.stateIdx.size(); i++) {

      vNoRecomb[i] = bSumState[hprior.stateIdx[i]];

    }

    for (size_t i = 0; i < hprior.nState(); i++) {

      tmpBw[i] = 0;

      for (size_t j = 0; j < lk.size(); j++) {

        tmpBw[i] += (lk[j] * bwd.back()[j]) * ibdRecombProbs.pRec_[siteI - 1];

      }

      tmpBw[i] *= statePrior[i];
      tmpBw[i] += lk[i] * (ibdRecombProbs.pNoRec_[siteI - 1]) * vNoRecomb[i];
      tmpBw[i] *= hprior.priorProb[i][siteI];

    }

    normalizeBySum(tmpBw);
    bwd.push_back(tmpBw);

  }

  reverse(bwd.begin(), bwd.end());

}


void kgd::IBDpath::computeIbdPathFwdProb(std::vector<double> proportion, std::vector<double> statePrior) {

  fm.clear();

  std::vector<double> vPrior = vecProd(statePrior, hprior.priorProbTrans[0]);

  std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, 0);

  updateFmAtSiteI(vPrior, lk);

  for (size_t siteI = 1; siteI < nLoci(); siteI++) {

    std::vector<double> vNoRec;

    for (size_t stateIdxTmp : hprior.stateIdx) {

      vNoRec.push_back(fSumState[stateIdxTmp]);

    }

    for (size_t i = 0; i < hprior.nState(); i++) {

      vPrior[i] = (vNoRec[i] * ibdRecombProbs.pNoRec_[siteI] + fSum * ibdRecombProbs.pRec_[siteI] * statePrior[i]) * hprior.priorProbTrans[siteI][i];

    }

    lk = computeLlkOfStatesAtSiteI(proportion, siteI);
    //cout << "lk = " ; for (double l :lk){cout << l << " ";}cout<<endl;

    updateFmAtSiteI(vPrior, lk);
    //for (double p : fm.back()){printf("%8.4f ", p);}cout<<endl;

  }

}


void kgd::IBDpath::updateFmAtSiteI(std::vector<double> &prior, std::vector<double> &llk) {

  std::vector<double> postAtSiteI = vecProd(prior, llk);
  //normalizeByMax(postAtSiteI);

  normalizeBySum(postAtSiteI);

  fm.push_back(postAtSiteI);
  fSum = sumOfVec(postAtSiteI);

  for (size_t i = 0; i < fSumState.size(); i++) {

    fSumState[i] = 0;

    for (size_t j = 0; j < hprior.nState(); j++) {

      fSumState[i] += ibdTransProbs[i][j] * postAtSiteI[j];

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

    std::vector<int> hSetI = hprior.hSet[indx];

    double qs = 0;

    for (size_t j = 0; j < kStrain(); j++) {

      qs += (double) hSetI[j] * proportion[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;

    if ((qs > 0) & (qs < 1)) {

      sumLLK += logBetaPdf(qs2, llkSurf[i][0], llkSurf[i][1]);

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

  assert(ibdTransProbs.size() == 0);

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

  return hprior.getIBDconfigureHeader();

}


std::vector<std::vector<double> > kgd::IBDpath::reshapeProbs(std::vector<std::vector<double> > &probs) {

  assert(nLoci() == probs.size());

  std::vector<std::vector<double> > ret;

  for (size_t siteIndex = 0; siteIndex < nLoci(); siteIndex++) {

    size_t previousStateIdx = 0;
    std::vector<double> tmpRow;
    double cumProb = 0;

    for (size_t prob_ij = 0; prob_ij < probs[siteIndex].size(); prob_ij++) {

      cumProb += probs[siteIndex][prob_ij];

      if (previousStateIdx != hprior.stateIdx[prob_ij]) {
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

  std::vector<double> pr0(kStrain());

  for (int i = 0; i < (int) pr0.size(); i++) {

    pr0[i] = binomialPdf(i, (int) (kStrain() - 1), theta);

  }

  std::vector<double> effectiveKPrior;

  for (size_t effectiveKtmp : hprior.effectiveK) {

    int effectiveKidx = effectiveKtmp - 1;
    assert(effectiveKidx >= 0);
    assert(effectiveKidx < (int) kStrain());
    effectiveKPrior.push_back(pr0[effectiveKidx] / uniqueEffectiveKCount[effectiveKidx]);

  }

  return effectiveKPrior;

}


std::vector<double> kgd::IBDpath::computeStatePrior(std::vector<double> effectiveKPrior) {

  std::vector<double> ret;

  for (size_t stateIdxTmp : hprior.stateIdx) {

    ret.push_back(effectiveKPrior[stateIdxTmp]);

  }

  return ret;

}


void kgd::IBDpath::computeAndUpdateTheta() {

  std::vector<size_t> obsState;

  size_t previousState = 0;
  size_t atSiteI = 0;

  for (size_t a : ibdConfigurePath) {

    if (a != previousState) {

      obsState.push_back(a);

    }

    if (hprior.stateIdx[a] != hprior.stateIdx[previousState]) {

      IBDpathChangeAt[atSiteI] += 1.0;
      currentIBDpathChangeAt[atSiteI] = 1.0;

    } else {

      currentIBDpathChangeAt[atSiteI] = 0.0;

    }

    previousState = a;
    atSiteI++;

  }

  size_t sumOfKeffStates = 0;
  size_t sccs = 0;

  for (size_t obs : obsState) {

    sumOfKeffStates += hprior.effectiveK[obs] - 1;
    sccs += kStrain() - hprior.effectiveK[obs];

  }
  //setTheta(rBeta(sccs+1.0, sumOfKeffStates+1.0, propRg_));

  setTheta(rBeta(sccs + 1.0, sumOfKeffStates + 1.0, ibdRg_));

}


void kgd::IBDpath::computeUniqueEffectiveKCount() {

  uniqueEffectiveKCount = std::vector<int>(kStrain());

  for (size_t effectiveKtmp : hprior.effectiveK) {

    int effectiveKidx = effectiveKtmp - 1;

    assert(effectiveKidx >= 0);

    uniqueEffectiveKCount[effectiveKidx]++;

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

  assert(llkSurf.size() == nLoci());

}


std::vector<double> kgd::IBDpath::computeLlkOfStatesAtSiteI(std::vector<double> proportion, size_t siteI, double err) {

  std::vector<double> llks;

  for (std::vector<int> hSetI : hprior.hSet) {

    double qs = 0;

    for (size_t j = 0; j < kStrain(); j++) {

      qs += (double) hSetI[j] * proportion[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;
    //cout << qs2 << endl;

    llks.push_back(logBetaPdf(qs2, llkSurf[siteI][0], llkSurf[siteI][1]));

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


