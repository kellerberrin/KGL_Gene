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

#include <random>
#include <stdio.h>
#include <limits>       // std::numeric_limits< double >::min()
#include <math.h>       // ceil
#include <kgl_exec_env.h>
#include "kgd_global.h"     // dout
#include "kgd_updateHap.h"
#include "kgd_mcmc.h"
#include "kgd_utility.h"

namespace kgl = kellerberrin::genome;
namespace kgd = kellerberrin::deploid;



// initialiseMCMCmachinery
kgd::McmcMachinery::McmcMachinery(std::shared_ptr<DEploidIO> dEploidIO,
                                  std::shared_ptr<McmcSample> mcmcSample,
                                  std::shared_ptr<RandomGenerator> randomGenerator,
                                  bool useIBD) {

  dEploidIO_ = dEploidIO;
  panel_ = dEploidIO->panel;
  mcmcSample_ = mcmcSample;
  seed_ = randomGenerator->seed();
  hapRg_ = randomGenerator;
  mcmcEventRg_ = randomGenerator;
  propRg_ = randomGenerator;
  initialHapRg_ = randomGenerator;

  if (useIBD == true) {

    calcMaxIteration(100, 10, 0.5);

  } else {

    calcMaxIteration(dEploidIO_->nMcmcSample_, dEploidIO_->mcmcMachineryRate_, dEploidIO_->mcmcBurn_);

  }

  MN_LOG_TITRE = 0.0;
  SD_LOG_TITRE = (useIBD == true) ? dEploidIO_->ibdSigma() : dEploidIO_->parameterSigma();
  PROP_SCALE = 40.0;

  stdNorm_ = new StandNormalRandomSample(seed_);

  setKstrain(dEploidIO_->kStrain());
  setNLoci(dEploidIO_->plaf_.size());
  initializeMcmcChain(useIBD);

}


kgd::McmcMachinery::~McmcMachinery() {

  if (stdNorm_) {

    delete stdNorm_;

  }

}


void kgd::McmcMachinery::calcMaxIteration(size_t nSample, size_t McmcMachineryRate, double burnIn) {

  burnIn_ = burnIn;
  McmcMachineryRate_ = McmcMachineryRate;
  maxIteration_ = (size_t) ceil((double) nSample * (double) McmcMachineryRate / (1.0 - burnIn_)) + 1;
  mcmcThresh_ = (size_t) ceil((double) nSample * (double) McmcMachineryRate * burnIn_ / (1.0 - burnIn_));

}


void kgd::McmcMachinery::initializeMcmcChain(bool useIBD) {
  // Initialization

  kgl::ExecEnv::log().info("############ initializeMcmcChain() ###########");

  initializeTitre();
  currentLogPriorTitre_ = calcLogPriorTitre(currentTitre_);
  initializeHap();
  initializeProp();
  initializeExpectedWsaf(); // This requires currentHap_ and currentProp_
  currentLLks_ = calcLLKs(dEploidIO_->refCount_,
                          dEploidIO_->altCount_,
                          currentExpectedWsaf_,
                          0,
                          currentExpectedWsaf_.size(),
                          dEploidIO_->scalingFactor());

  acceptUpdate = 0;

  if (dEploidIO_->doAllowInbreeding() == true) {

    initializeUpdateReferencePanel(panel_->truePanelSize() + kStrain_ - 1);

  }

  if (useIBD == true) {

    ibdInitializeEssentials();

  }

  mcmcSample_->siteOfTwoSwitchOne = std::vector<double>(nLoci());
  mcmcSample_->siteOfTwoMissCopyOne = std::vector<double>(nLoci());
  mcmcSample_->siteOfTwoSwitchTwo = std::vector<double>(nLoci());
  mcmcSample_->siteOfTwoMissCopyTwo = std::vector<double>(nLoci());
  mcmcSample_->siteOfOneSwitchOne = std::vector<double>(nLoci());
  mcmcSample_->siteOfOneMissCopyOne = std::vector<double>(nLoci());

  mcmcSample_->currentsiteOfTwoSwitchOne = std::vector<double>(nLoci());
  mcmcSample_->currentsiteOfTwoMissCopyOne = std::vector<double>(nLoci());
  mcmcSample_->currentsiteOfTwoSwitchTwo = std::vector<double>(nLoci());
  mcmcSample_->currentsiteOfTwoMissCopyTwo = std::vector<double>(nLoci());
  mcmcSample_->currentsiteOfOneSwitchOne = std::vector<double>(nLoci());
  mcmcSample_->currentsiteOfOneMissCopyOne = std::vector<double>(nLoci());

  assert (doutProp());
  assert (doutLLK());

  kgl::ExecEnv::log().info("############ initializeMcmcChain() Complete ###########");

}


void kgd::McmcMachinery::initializeHap() {

  assert(currentHap_.size() == 0);

  if (dEploidIO_->initialHapWasGiven()) {

    currentHap_ = dEploidIO_->initialHap_;

  } else {

    for (size_t i = 0; i < dEploidIO_->plaf_.size(); i++) {

      double currentPlaf = dEploidIO_->plaf_[i];
      std::vector<double> tmpVec;

      for (size_t k = 0; k < kStrain_; k++) {

        tmpVec.push_back(rBernoulli(currentPlaf));

      }

      currentHap_.push_back(tmpVec);

    }

  }

  assert(currentHap_.size() == dEploidIO_->plaf_.size());

}


void kgd::McmcMachinery::initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo) {

  if (dEploidIO_->doAllowInbreeding() != true) {

    return;

  }

  panel_->initializeUpdatePanel(inbreedingPanelSizeSetTo);

}


void kgd::McmcMachinery::updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain) {

  if (burnIn_ > currentMcmcIteration_) {

    return;

  }

  //if ( dEploidIO_->doAllowInbreeding() != true ){
  //return;
  //}
  panel_->updatePanelWithHaps(inbreedingPanelSizeSetTo, excludedStrain, currentHap_);
}


double kgd::McmcMachinery::rBernoulli(double p) {

  double u = initialHapRg_->sample();
  return (u < p) ? 1.0 : 0.0;

}


void kgd::McmcMachinery::initializeExpectedWsaf() {

  assert(currentExpectedWsaf_.size() == 0);

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

  assert(currentExpectedWsaf_.size() == nLoci_);

  cumExpectedWsaf_ = currentExpectedWsaf_;

}


void kgd::McmcMachinery::initializellk() {

  assert(currentLLks_.size() == (size_t) 0);

  currentLLks_ = std::vector<double>(nLoci_, 0.0);

  assert(currentLLks_.size() == nLoci_);

}


void kgd::McmcMachinery::initializeProp() {

  assert(currentProp_.size() == (size_t) 0);

  currentProp_ = (dEploidIO_->initialPropWasGiven()) ? dEploidIO_->initialProp_ : titre2prop(currentTitre_);

  if (dEploidIO_->initialPropWasGiven()) {

    currentTitre_.clear();

    for (size_t i = 0; i < dEploidIO_->initialProp_.size(); i++) {

      currentTitre_.push_back(log(dEploidIO_->initialProp_[i]));

    }

  }

}


double kgd::McmcMachinery::calcLogPriorTitre(std::vector<double> &tmpTitre) {
  //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));

  std::vector<double> tmp;

  for (auto const &value: tmpTitre) {

    tmp.push_back(log(normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)));

  }

  return sumOfVec(tmp);

}


void kgd::McmcMachinery::initializeTitre() {
  /*   titre<-rnorm(initial.k, MN_LOG_TITRE, SD_LOG_TITRE); */

  assert(currentTitre_.size() == 0);

  currentTitre_ = std::vector<double>(kStrain_, 0.0);

  if (dEploidIO_->doUpdateProp()) {

    for (size_t k = 0; k < kStrain_; k++) {

      currentTitre_[k] = initialTitreNormalVariable();

    }

  }

  assert(currentTitre_.size() == kStrain_);

}


std::vector<double> kgd::McmcMachinery::titre2prop(std::vector<double> &tmpTitre) {

  std::vector<double> tmpExpTitre;

  for (auto const &value: tmpTitre) {

    tmpExpTitre.push_back(exp(value));

  }

  double tmpSum = sumOfVec(tmpExpTitre);

  std::vector<double> tmpProp;

  for (auto const &value: tmpExpTitre) {

    tmpProp.push_back(value / tmpSum);
    assert (tmpProp.back() > 0);
    assert (tmpProp.back() <= 1);

  }

  return tmpProp;

}


void kgd::McmcMachinery::initializePropIBD() {

  currentProp_ = (dEploidIO_->initialPropWasGiven()) ? dEploidIO_->initialProp_ : titre2prop(currentTitre_);

}


void kgd::McmcMachinery::runMcmcChain(bool showProgress, bool useIBD, bool notInR) {

  for (currentMcmcIteration_ = 0; currentMcmcIteration_ < maxIteration_; currentMcmcIteration_++) {

    kgl::ExecEnv::log().info("MCMC iteration: {}/{}", currentMcmcIteration_, maxIteration_);

    if (currentMcmcIteration_ > 0 && currentMcmcIteration_ % 100 == 0 && showProgress) {

#ifndef RBUILD
      std::clog << "\r" << " MCMC step" << std::setw(4) << int(currentMcmcIteration_ * 100 / maxIteration_) << "% completed." << std::flush;
#endif

    }

    sampleMcmcEvent(useIBD);
  }

#ifndef RBUILD
  std::clog << "\r" << " MCMC step" << std::setw(4) << 100 << "% completed." << std::endl;
#endif

  mcmcSample_->hap = currentHap_;

  writeLastFwdProb(useIBD);

  dEploidIO_->finalProp_ = mcmcSample_->proportion.back();

  for (size_t atSiteI = 0; atSiteI < nLoci(); atSiteI++) {

    mcmcSample_->siteOfTwoSwitchOne[atSiteI] /= (double) maxIteration_;
    mcmcSample_->siteOfTwoMissCopyOne[atSiteI] /= (double) maxIteration_;
    mcmcSample_->siteOfTwoSwitchTwo[atSiteI] /= (double) maxIteration_;
    mcmcSample_->siteOfTwoMissCopyTwo[atSiteI] /= (double) maxIteration_;
    mcmcSample_->siteOfOneSwitchOne[atSiteI] /= (double) maxIteration_;
    mcmcSample_->siteOfOneMissCopyOne[atSiteI] /= (double) maxIteration_;

  }

  if (notInR) {

    dEploidIO_->writeMcmcRelated(mcmcSample_, useIBD);

  }

  if (useIBD == true) {

    for (size_t atSiteI = 0; atSiteI < nLoci(); atSiteI++) {

      ibdPath.IBDpathChangeAt[atSiteI] /= (double) maxIteration_;

    }

    std::clog << "Proportion update acceptance rate: " << acceptUpdate / (kStrain() * 1.0 * maxIteration_) << std::endl;

    dEploidIO_->initialProp_ = averageProportion(mcmcSample_->proportion);
    dEploidIO_->setInitialPropWasGiven(true);
    dEploidIO_->setDoUpdateProp(false);
    dEploidIO_->initialHap_ = mcmcSample_->hap;
    dEploidIO_->setInitialHapWasGiven(true);

  }

  computeDiagnostics();


  dout << "###########################################" << std::endl;
  dout << "#            MCMC RUN finished            #" << std::endl;
  dout << "###########################################" << std::endl;
}


void kgd::McmcMachinery::computeDiagnostics() {

  //clog << "Proportion update acceptance rate: "<<acceptUpdate / (kStrain()*1.0*maxIteration_)<<endl;

  dEploidIO_->setacceptRatio(acceptUpdate / (1.0 * maxIteration_));

  // average cumulate expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); i++) {

    cumExpectedWsaf_[i] /= dEploidIO_->nMcmcSample_;

  }

  std::vector<double> tmpLLKs1 = calcLLKs(dEploidIO_->refCount_,
                                          dEploidIO_->altCount_,
                                          cumExpectedWsaf_,
                                          0,
                                          cumExpectedWsaf_.size(),
                                          dEploidIO_->scalingFactor());

  dEploidIO_->setmeanThetallks(sumOfVec(tmpLLKs1));

  std::vector<double> wsaf_vec;

  for (size_t i = 0; i < nLoci(); i++) {

    double wsaf = dEploidIO_->altCount_[i] / (dEploidIO_->refCount_[i] + dEploidIO_->altCount_[i] + 0.00000000000001);
    double adjustedWsaf = wsaf * (1 - 0.01) + (1 - wsaf) * 0.01;
    wsaf_vec.push_back(adjustedWsaf);
    //llkOfData.push_back( logBetaPdf(adjustedWsaf, llkSurf[i][0], llkSurf[i][1]));

  }

  std::vector<double> tmpLLKs = calcLLKs(dEploidIO_->refCount_,
                                         dEploidIO_->altCount_,
                                         wsaf_vec,
                                         0,
                                         wsaf_vec.size(),
                                         dEploidIO_->scalingFactor());

  dEploidIO_->setmaxLLKs(sumOfVec(tmpLLKs));

  double sum = std::accumulate(mcmcSample_->sumLLKs.begin(), mcmcSample_->sumLLKs.end(), 0.0);
  double mean = sum / mcmcSample_->sumLLKs.size();
  double sq_sum = std::inner_product(mcmcSample_->sumLLKs.begin(), mcmcSample_->sumLLKs.end(),
                                     mcmcSample_->sumLLKs.begin(), 0.0);
  double varLLKs = sq_sum / mcmcSample_->sumLLKs.size() - mean * mean;
  double stdev = std::sqrt(varLLKs);

  dEploidIO_->setmeanllks(mean);
  dEploidIO_->setstdvllks(stdev);

  double dicByVar = (-2 * mean) + 4 * varLLKs / 2;
  dEploidIO_->setdicByVar(dicByVar);
  //return (  mean(-2*tmpllk) + var(-2*tmpllk)/2 )# D_bar + 1/2 var (D_theta), where D_theta = -2*tmpllk, and D_bar = mean(D_theta)

  double dicWSAFBar = -2 * sumOfVec(tmpLLKs1);
  double dicByTheta = (-2 * mean) + (-2 * mean) - dicWSAFBar;
  dEploidIO_->setdicByTheta(dicByTheta);
  //DIC.WSAF.bar = -2 * sum(thetallk)
  //return (  mean(-2*tmpllk) + (mean(-2*tmpllk) - DIC.WSAF.bar) ) # D_bar + pD, where pD = D_bar - D_theta, and D_bar = mean(D_theta)

}

std::vector<double> kgd::McmcMachinery::averageProportion(std::vector<std::vector<double> > &proportion) {

  assert(proportion.size() > 0);

  std::vector<double> ret(kStrain());

  for (size_t i = 0; i < kStrain(); i++) {

    for (std::vector<double> p : proportion) {

      ret[i] += p[i];

    }

    ret[i] /= (1.0 * proportion.size());

  }

  normalizeBySum(ret);

  return ret;
}


void kgd::McmcMachinery::sampleMcmcEvent(bool useIBD) {

  recordingMcmcBool_ = (currentMcmcIteration_ > mcmcThresh_ && currentMcmcIteration_ % McmcMachineryRate_ == 0);

  if (useIBD == true) {

    ibdSampleMcmcEventStep();
    assert(doutProp());

  } else {

    eventInt_ = mcmcEventRg_->sampleInt(3);

    if ((eventInt_ == 0) && (dEploidIO_->doUpdateProp())) {
      updateProportion();

    } else if ((eventInt_ == 1) && (dEploidIO_->doUpdateSingle())) {

      updateSingleHap();

    } else if ((eventInt_ == 2) && (dEploidIO_->doUpdatePair())) {

      updatePairHaps();

    }

  }

  assert(doutLLK());

  if (recordingMcmcBool_) {
    recordMcmcMachinery();
  }

}


void kgd::McmcMachinery::ibdInitializeEssentials() {

  initializePropIBD();
  ibdPath.init(*dEploidIO_, hapRg_);

  std::vector<double> llkOfData;

  for (size_t i = 0; i < nLoci(); i++) {

    double wsaf = dEploidIO_->altCount_[i] / (dEploidIO_->refCount_[i] + dEploidIO_->altCount_[i] + 0.00000000000001);
    double adjustedWsaf = wsaf * (1 - 0.01) + (1 - wsaf) * 0.01;
    llkOfData.push_back(logBetaPdf(adjustedWsaf, ibdPath.llkSurf[i][0], ibdPath.llkSurf[i][1]));

  }

  dout << "LLK of data = " << sumOfVec(llkOfData) << std::endl;

}


void kgd::McmcMachinery::ibdSampleMcmcEventStep() {

  std::vector<double> effectiveKPrior = ibdPath.computeEffectiveKPrior(ibdPath.theta());
  std::vector<double> statePrior = ibdPath.computeStatePrior(effectiveKPrior);
  // First building the path likelihood
  ibdPath.computeIbdPathFwdProb(currentProp_, statePrior);

  ////#Now sample path given matrix
  ibdPath.ibdSamplePath(statePrior);

  //#Get haplotypes and update LLK for each site
  ibdUpdateHaplotypesFromPrior();
  std::vector<double> llkAtAllSites = computeLlkAtAllSites();
  ////#Given current haplotypes, sample titres 1 by 1 using MH
  ibdUpdateProportionGivenHap(llkAtAllSites);
  // Compute new theta after all proportion and haplotypes are up to date.
  ibdPath.computeAndUpdateTheta();

  currentLLks_ = llkAtAllSites;
  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::McmcMachinery::ibdUpdateHaplotypesFromPrior() {

  for (size_t i = 0; i < nLoci(); i++) {

    for (size_t j = 0; j < kStrain(); j++) {

      currentHap_[i][j] = (double) ibdPath.hprior.hSet[ibdPath.ibdConfigurePath[i]][j];

    }

  }

}


void kgd::McmcMachinery::ibdUpdateProportionGivenHap(std::vector<double> &llkAtAllSites) {

  for (size_t i = 0; i < kStrain(); i++) {

    double v0 = currentTitre_[i];

    std::vector<double> oldProp = currentProp_;
    //currentTitre_[i] += (stdNorm_->genReal() * 0.1 + 0.0); // tit.0[i]+rnorm(1, 0, scale.t.prop);

    currentTitre_[i] += (stdNorm_->genReal() * SD_LOG_TITRE * 1.0 / PROP_SCALE + 0.0); // tit.0[i]+rnorm(1, 0, scale.t.prop);
    currentProp_ = titre2prop(currentTitre_);
    std::vector<double> vv = computeLlkAtAllSites();

    double rr = normal_pdf(currentTitre_[i], 0, 1) /
                normal_pdf(v0, 0, 1) * exp(sumOfVec(vv) - sumOfVec(llkAtAllSites));

    if (propRg_->sample() < rr) {

      llkAtAllSites = vv;
      acceptUpdate++;

    } else {

      currentTitre_[i] = v0;
      currentProp_ = oldProp;

    }

  }

}


std::vector<double> kgd::McmcMachinery::computeLlkAtAllSites(double err) {

  std::vector<double> ret;

  for (size_t site = 0; site < nLoci(); site++) {

    double qs = 0;

    for (size_t j = 0; j < kStrain(); j++) {

      qs += (double) currentHap_[site][j] * currentProp_[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;

    ret.push_back(logBetaPdf(qs2, ibdPath.llkSurf[site][0], ibdPath.llkSurf[site][1]));

  }

  return ret;

}


std::vector<double> kgd::McmcMachinery::calcExpectedWsaf(std::vector<double> &proportion) {
  //assert ( sumOfVec(proportion) == 1.0); // this fails ...

  std::vector<double> expectedWsaf(nLoci_, 0.0);

  for (size_t i = 0; i < currentHap_.size(); i++) {

    assert(kStrain_ == currentHap_[i].size());

    for (size_t k = 0; k < kStrain_; k++) {

      expectedWsaf[i] += currentHap_[i][k] * proportion[k];

    }

    assert (expectedWsaf[i] >= 0);
    //assert ( expectedWsaf[i] <= 1.0 );
  }

  return expectedWsaf;

}


void kgd::McmcMachinery::recordMcmcMachinery() {

  dout << "***Record mcmc sample " << std::endl;

  mcmcSample_->proportion.push_back(currentProp_);
  mcmcSample_->sumLLKs.push_back(sumOfVec(currentLLks_));
  mcmcSample_->moves.push_back(eventInt_);

  // Cumulate expectedWSAF for computing the mean expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); i++) {

    cumExpectedWsaf_[i] += currentExpectedWsaf_[i];

  }

}


void kgd::McmcMachinery::updateProportion() {

  dout << " Attempt of update proportion";

  if (kStrain_ < 2) {

    dout << "(failed)" << std::endl;
    return;

  }

  // calculate dt
  std::vector<double> tmpTitre = calcTmpTitre();
  std::vector<double> tmpProp = titre2prop(tmpTitre);

  if (min_value(tmpProp) < 0 || max_value(tmpProp) > 1) {

    dout << "(failed)" << std::endl;
    return;

  }

  std::vector<double> tmpExpecedWsaf = calcExpectedWsaf(tmpProp);
  std::vector<double> tmpLLKs = calcLLKs(dEploidIO_->refCount_,
                                         dEploidIO_->altCount_,
                                         tmpExpecedWsaf, 0,
                                         tmpExpecedWsaf.size(),
                                         dEploidIO_->scalingFactor());

  double diffLLKs = deltaLLKs(tmpLLKs);
  double tmpLogPriorTitre = calcLogPriorTitre(tmpTitre);
  double priorPropRatio = exp(tmpLogPriorTitre - currentLogPriorTitre_);
  double hastingsRatio = 1.0;

  //runif(1)<prior.prop.ratio*hastings.ratio*exp(del.llk))
  if (propRg_->sample() > priorPropRatio * hastingsRatio * exp(diffLLKs)) {

    dout << "(failed)" << std::endl;
    return;

  }

  dout << "(successed) " << std::endl;
  acceptUpdate++;

  currentExpectedWsaf_ = tmpExpecedWsaf;
  currentLLks_ = tmpLLKs;
  currentLogPriorTitre_ = tmpLogPriorTitre;
  currentTitre_ = tmpTitre;
  currentProp_ = tmpProp;

  assert (doutProp());

}


double kgd::McmcMachinery::deltaLLKs(std::vector<double> &newLLKs) {

  std::vector<double> tmpdiff = vecDiff(newLLKs, currentLLks_);

  return sumOfVec(tmpdiff);

}


std::vector<double> kgd::McmcMachinery::calcTmpTitre() {

  std::vector<double> tmpTitre;

  for (size_t k = 0; k < kStrain_; k++) {

    double dt = deltaXnormalVariable();
    tmpTitre.push_back(currentTitre_[k] + dt);

  }

  return tmpTitre;

}


void kgd::McmcMachinery::updateSingleHap() {

  kgl::ExecEnv::log().info("McmcMachinery::updateSingleHap()");

  findUpdatingStrainSingle();

  if (dEploidIO_->doAllowInbreeding() == true) {

    updateReferencePanel(panel_->truePanelSize() + kStrain_ - 1, strainIndex_);

  }


  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts_.size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts_[chromi];
    size_t length = dEploidIO_->position_[chromi].size();

    kgl::ExecEnv::log().info("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdateSingleHap updating(dEploidIO_->refCount_,
                             dEploidIO_->altCount_,
                             dEploidIO_->plaf_,
                             currentExpectedWsaf_,
                             currentProp_,
                             currentHap_,
                             hapRg_,
                             start,
                             length,
                             panel_,
                             dEploidIO_->missCopyProb_,
                             dEploidIO_->scalingFactor(),
                             strainIndex_);

    if (dEploidIO_->doAllowInbreeding()) {

      updating.setPanelSize(panel_->inbreedingPanelSize());

    }

    updating.core(dEploidIO_->refCount_,
                  dEploidIO_->altCount_,
                  dEploidIO_->plaf_,
                  currentExpectedWsaf_,
                  currentProp_,
                  currentHap_);

    size_t updateIndex = 0;

    for (size_t ii = start; ii < (start + length); ii++) {

      currentHap_[ii][strainIndex_] = updating.hap_[updateIndex];
      currentLLks_[ii] = updating.newLLK[updateIndex];
      updateIndex++;

    }

    for (size_t siteI = 0; siteI < length; siteI++) {

      mcmcSample_->siteOfOneSwitchOne[start + siteI] += updating.siteOfOneSwitchOne[siteI];
      mcmcSample_->siteOfOneMissCopyOne[start + siteI] += updating.siteOfOneMissCopyOne[siteI];
      mcmcSample_->siteOfOneSwitchOne[start + siteI] = updating.siteOfOneSwitchOne[siteI];
      mcmcSample_->siteOfOneMissCopyOne[start + siteI] = updating.siteOfOneMissCopyOne[siteI];

    }

  }

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::McmcMachinery::updatePairHaps() {

  if (kStrain() == 1) {
    return;
  }

  dout << " Update Pair Hap " << std::endl;
  findUpdatingStrainPair();

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts_.size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts_[chromi];
    size_t length = dEploidIO_->position_[chromi].size();

    dout << "   Update Chrom with index " << chromi << ", starts at " << start << ", with " << length << " sites" << std::endl;

    UpdatePairHap updating(dEploidIO_->refCount_,
                           dEploidIO_->altCount_,
                           dEploidIO_->plaf_,
                           currentExpectedWsaf_,
                           currentProp_,
                           currentHap_,
                           hapRg_,
                           start,
                           length,
                           panel_,
                           dEploidIO_->missCopyProb_,
                           dEploidIO_->scalingFactor(),
                           dEploidIO_->forbidCopyFromSame(),
                           strainIndex1_,
                           strainIndex2_);

    updating.core(dEploidIO_->refCount_,
                  dEploidIO_->altCount_,
                  dEploidIO_->plaf_,
                  currentExpectedWsaf_,
                  currentProp_,
                  currentHap_);

    size_t updateIndex = 0;
    for (size_t ii = start; ii < (start + length); ii++) {

      currentHap_[ii][strainIndex1_] = updating.hap1_[updateIndex];
      currentHap_[ii][strainIndex2_] = updating.hap2_[updateIndex];
      currentLLks_[ii] = updating.newLLK[updateIndex];
      updateIndex++;

    }

    for (size_t siteI = 0; siteI < length; siteI++) {

      mcmcSample_->siteOfTwoSwitchOne[start + siteI] += updating.siteOfTwoSwitchOne[siteI];
      mcmcSample_->siteOfTwoMissCopyOne[start + siteI] += updating.siteOfTwoMissCopyOne[siteI];
      mcmcSample_->siteOfTwoSwitchTwo[start + siteI] += updating.siteOfTwoSwitchTwo[siteI];
      mcmcSample_->siteOfTwoMissCopyTwo[start + siteI] += updating.siteOfTwoMissCopyTwo[siteI];
      mcmcSample_->currentsiteOfTwoSwitchOne[start + siteI] = updating.siteOfTwoSwitchOne[siteI];
      mcmcSample_->currentsiteOfTwoMissCopyOne[start + siteI] = updating.siteOfTwoMissCopyOne[siteI];
      mcmcSample_->currentsiteOfTwoSwitchTwo[start + siteI] = updating.siteOfTwoSwitchTwo[siteI];
      mcmcSample_->currentsiteOfTwoMissCopyTwo[start + siteI] = updating.siteOfTwoMissCopyTwo[siteI];

    }
  }

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::McmcMachinery::findUpdatingStrainSingle() {

  std::vector<double> eventProb(kStrain_, 1);

  normalizeBySum(eventProb);
  strainIndex_ = sampleIndexGivenProp(mcmcEventRg_, eventProb);

}


void kgd::McmcMachinery::findUpdatingStrainPair() {

  std::vector<size_t> strainIndex(2, 0);

  int t = 0; // total input records dealt with
  int m = 0; // number of items selected so far
  double u;

  while (m < 2) {

    u = mcmcEventRg_->sample(); // call a uniform(0,1) kgd_random number generator

    if ((kStrain_ - t) * u < 2 - m) {
      strainIndex[m] = t;
      m++;
    }

    t++;

  }

  strainIndex1_ = strainIndex[0];
  strainIndex2_ = strainIndex[1];

  assert(strainIndex1_ != strainIndex2_);
  dout << "  Updating hap: " << strainIndex1_ << " and " << strainIndex2_ << std::endl;

}


