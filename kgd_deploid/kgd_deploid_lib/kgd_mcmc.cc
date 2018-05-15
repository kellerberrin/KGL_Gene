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
#include "kgd_update_haplotype.h"
#include "kgd_update_single_haplotype.h"  // chromPainting
#include "kgd_update_pair_haplotype.h"
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
  panel_ = dEploidIO->getPanel();
  mcmcSample_ = mcmcSample;
  seed_ = randomGenerator->seed();
  hapRg_ = randomGenerator;
  mcmcEventRg_ = randomGenerator;
  propRg_ = randomGenerator;
  initialHapRg_ = randomGenerator;

  if (useIBD == true) {

    calcMaxIteration(100, 10, 0.5);

  } else {

    calcMaxIteration(dEploidIO_->getMcmcSample(), dEploidIO_->getMcmcMachineryRate(), dEploidIO_->getMcmcBurn());

  }

  MN_LOG_TITRE = 0.0;
  SD_LOG_TITRE = (useIBD == true) ? dEploidIO_->ibdSigma() : dEploidIO_->parameterSigma();
  PROP_SCALE = 40.0;

  stdNorm_ = std::make_shared<StandNormalRandomSample>(seed_);

  setKstrain(dEploidIO_->kStrain());
  setNLoci(dEploidIO_->getPlaf().size());
  initializeMcmcChain(useIBD);

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
  currentLLks_ = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                   dEploidIO_->getAltCount(),
                                   currentExpectedWsaf_,
                                   0,
                                   currentExpectedWsaf_.size(),
                                   dEploidIO_->scalingFactor());

  acceptUpdate = 0;

  if (dEploidIO_->doAllowInbreeding()) {

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

}


void kgd::McmcMachinery::initializeHap() {

  assert(currentHap_.size() == 0);

  if (dEploidIO_->initialHapWasGiven()) {

    currentHap_ = dEploidIO_->getInitialHap();

  } else {

    for (size_t i = 0; i < dEploidIO_->getPlaf().size(); i++) {

      double currentPlaf = dEploidIO_->getPlaf()[i];
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

  if (dEploidIO_->doAllowInbreeding()) {

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

  currentProp_ = (dEploidIO_->initialPropWasGiven()) ? dEploidIO_->getInitialProp() : titre2prop(currentTitre_);

  if (dEploidIO_->initialPropWasGiven()) {

    currentTitre_.clear();

    for (size_t i = 0; i < dEploidIO_->getInitialProp().size(); i++) {

      currentTitre_.push_back(log(dEploidIO_->getInitialProp()[i]));

    }

  }

}


double kgd::McmcMachinery::calcLogPriorTitre(std::vector<double> &tmpTitre) {
  //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));

  std::vector<double> tmp;

  for (auto const &value: tmpTitre) {

    tmp.push_back(log(Utility::normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)));

  }

  return Utility::sumOfVec(tmp);

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

  double tmpSum = Utility::sumOfVec(tmpExpTitre);

  std::vector<double> tmpProp;

  for (auto const &value: tmpExpTitre) {

    tmpProp.push_back(value / tmpSum);
    assert (tmpProp.back() > 0);
    assert (tmpProp.back() <= 1);

  }

  return tmpProp;

}


void kgd::McmcMachinery::initializePropIBD() {

  currentProp_ = (dEploidIO_->initialPropWasGiven()) ? dEploidIO_->getInitialProp() : titre2prop(currentTitre_);

}


void kgd::McmcMachinery::runMcmcChain(bool showProgress, bool useIBD, bool notInR) {

  for (currentMcmcIteration_ = 0; currentMcmcIteration_ < maxIteration_; currentMcmcIteration_++) {

    if (currentMcmcIteration_ > 0 && currentMcmcIteration_ % 100 == 0 && showProgress) {

      kgl::ExecEnv::log().info("MCMC iteration: {}/{}, completed: {}%", currentMcmcIteration_, maxIteration_, int(currentMcmcIteration_ * 100 / maxIteration_));

    }

    sampleMcmcEvent(useIBD);

  }

  mcmcSample_->hap = currentHap_;

  writeLastFwdProb(useIBD);

  dEploidIO_->setFinalProp(mcmcSample_->proportion.back());

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

      ibdPath.IBDPathChangeAt(atSiteI,  static_cast<double>(maxIteration_));

    }

    kgl::ExecEnv::log().info("Proportion update acceptance rate: {}", acceptUpdate / (kStrain() * 1.0 * maxIteration_));

    dEploidIO_->setInitialProp(averageProportion(mcmcSample_->proportion));
    dEploidIO_->setInitialPropWasGiven(true);
    dEploidIO_->setDoUpdateProp(false);
    dEploidIO_->setInitialHap(mcmcSample_->hap);
    dEploidIO_->setInitialHapWasGiven(true);

  }

  computeDiagnostics();

  kgl::ExecEnv::log().info("#### MCMC RUN finished ####");

}


void kgd::McmcMachinery::computeDiagnostics() {

  dEploidIO_->setacceptRatio(acceptUpdate / (1.0 * maxIteration_));

  // average cumulate expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); i++) {

    cumExpectedWsaf_[i] /= dEploidIO_->getMcmcSample();

  }

  std::vector<double> tmpLLKs1 = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                                   dEploidIO_->getAltCount(),
                                                   cumExpectedWsaf_,
                                                   0,
                                                   cumExpectedWsaf_.size(),
                                                   dEploidIO_->scalingFactor());

  dEploidIO_->setmeanThetallks(Utility::sumOfVec(tmpLLKs1));

  std::vector<double> wsaf_vec;

  for (size_t i = 0; i < nLoci(); i++) {

    double wsaf = dEploidIO_->getAltCount()[i] / (dEploidIO_->getRefCount()[i] + dEploidIO_->getAltCount()[i] + 0.00000000000001);
    double adjustedWsaf = wsaf * (1 - 0.01) + (1 - wsaf) * 0.01;

    wsaf_vec.push_back(adjustedWsaf);

    //llkOfData.push_back( logBetaPdf(adjustedWsaf, llk_surf_[i][0], llk_surf_[i][1]));

  }

  std::vector<double> tmpLLKs = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                                  dEploidIO_->getAltCount(),
                                                  wsaf_vec,
                                                  0,
                                                  wsaf_vec.size(),
                                                  dEploidIO_->scalingFactor());

  dEploidIO_->setmaxLLKs(Utility::sumOfVec(tmpLLKs));

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

  double dicWSAFBar = -2 * Utility::sumOfVec(tmpLLKs1);
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

  Utility::normalizeBySum(ret);

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

    double wsaf = dEploidIO_->getAltCount()[i] / (dEploidIO_->getRefCount()[i] + dEploidIO_->getAltCount()[i] + 0.00000000000001);

    double adjustedWsaf = wsaf * (1 - 0.01) + (1 - wsaf) * 0.01;

    llkOfData.push_back(Utility::logBetaPdf(adjustedWsaf, ibdPath.getLogLikelihoodSurface()[i][0], ibdPath.getLogLikelihoodSurface()[i][1]));

  }

  dout << "LLK of data = " << Utility::sumOfVec(llkOfData) << std::endl;

}


void kgd::McmcMachinery::ibdSampleMcmcEventStep() {

  // Update the idb path.
  ibdPath.McmcUpdateStep(currentProp_);

  //#Get haplotypes and update LLK for each site
  ibdUpdateHaplotypesFromPrior();

  std::vector<double> llkAtAllSites = computeLlkAtAllSites();
  ////#Given current haplotypes, sample titres 1 by 1 using MH

  ibdUpdateProportionGivenHap(llkAtAllSites);

  currentLLks_ = llkAtAllSites;

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::McmcMachinery::ibdUpdateHaplotypesFromPrior() {

  for (size_t loci = 0; loci < nLoci(); ++loci) {

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      currentHap_[loci][strain] = ibdPath.UpdateHaplotypesFromPrior(strain, loci);

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

    double rr = Utility::normal_pdf(currentTitre_[i], 0, 1) /
    Utility::normal_pdf(v0, 0, 1) * exp(Utility::sumOfVec(vv) - Utility::sumOfVec(llkAtAllSites));

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

  for (size_t loci = 0; loci < nLoci(); ++loci) {

    double qs = 0;

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      qs += (double) currentHap_[loci][strain] * currentProp_[strain];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;

    ret.push_back(Utility::logBetaPdf(qs2, ibdPath.getLogLikelihoodSurface()[loci][0], ibdPath.getLogLikelihoodSurface()[loci][1]));

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
  mcmcSample_->sumLLKs.push_back(Utility::sumOfVec(currentLLks_));
  mcmcSample_->moves.push_back(eventInt_);

  // Cumulate expectedWSAF for computing the mean expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); i++) {

    cumExpectedWsaf_[i] += currentExpectedWsaf_[i];

  }

}


void kgd::McmcMachinery::updateProportion() {

  kgl::ExecEnv::log().vinfo("MCMC Attempt of update proportion");

  if (kStrain_ < 2) {

    kgl::ExecEnv::log().vinfo("MCMC Update failed (only 1 strain)");
    return;

  }

  // calculate dt
  std::vector<double> tmpTitre = calcTmpTitre();
  std::vector<double> tmpProp = titre2prop(tmpTitre);

  if (Utility::min_value(tmpProp) < 0 || Utility::max_value(tmpProp) > 1) {

    kgl::ExecEnv::log().vinfo("MCMC Update failed, tmpProp < 0 or tmpProp > 1");
    return;

  }

  std::vector<double> tmpExpecedWsaf = calcExpectedWsaf(tmpProp);

  std::vector<double> tmpLLKs = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                                  dEploidIO_->getAltCount(),
                                                  tmpExpecedWsaf,
                                                  0,
                                                  tmpExpecedWsaf.size(),
                                                  dEploidIO_->scalingFactor());

  double diffLLKs = deltaLLKs(tmpLLKs);

  double tmpLogPriorTitre = calcLogPriorTitre(tmpTitre);

  double priorPropRatio = exp(tmpLogPriorTitre - currentLogPriorTitre_);

  double hastingsRatio = 1.0;

  double sample_value = propRg_->sample();

  double priorRatio = priorPropRatio * hastingsRatio * exp(diffLLKs);

  if (sample_value > priorRatio) {

    kgl::ExecEnv::log().vinfo("MCMC Update failed, sample: {}, priorRatio :{}", sample_value, priorRatio);
    return;

  }

  kgl::ExecEnv::log().vinfo("MCMC Update succeeded, sample: {}, priorRatio :{}", sample_value, priorRatio);

  acceptUpdate++;
  currentExpectedWsaf_ = tmpExpecedWsaf;
  currentLLks_ = tmpLLKs;
  currentLogPriorTitre_ = tmpLogPriorTitre;
  currentTitre_ = tmpTitre;
  currentProp_ = tmpProp;

  assert (doutProp());

}


double kgd::McmcMachinery::deltaLLKs(std::vector<double> &newLLKs) {

  std::vector<double> tmpdiff = Utility::vecDiff(newLLKs, currentLLks_);

  return Utility::sumOfVec(tmpdiff);

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

  kgl::ExecEnv::log().vinfo("McmcMachinery::updateSingleHap()");

  findUpdatingStrainSingle();

  if (dEploidIO_->doAllowInbreeding()) {

    updateReferencePanel(panel_->truePanelSize() + kStrain_ - 1, strainIndex_);

  }

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts()[chromi];
    size_t length = dEploidIO_->getIndexPosition(chromi).size();

    kgl::ExecEnv::log().vinfo("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdateSingleHap updating(dEploidIO_->getRefCount(),
                             dEploidIO_->getAltCount(),
                             dEploidIO_->getPlaf(),
                             currentExpectedWsaf_,
                             currentProp_,
                             currentHap_,
                             hapRg_,
                             start,
                             length,
                             panel_,
                             dEploidIO_->getMissCopyProb(),
                             dEploidIO_->scalingFactor(),
                             strainIndex_);

    if (dEploidIO_->doAllowInbreeding()) {

      updating.setPanelSize(panel_->inbreedingPanelSize());

    }

    updating.core(dEploidIO_->getRefCount(),
                  dEploidIO_->getAltCount(),
                  dEploidIO_->getPlaf(),
                  currentExpectedWsaf_,
                  currentProp_,
                  currentHap_);

    size_t updateIndex = 0;

    for (size_t ii = start; ii < (start + length); ii++) {

      currentHap_[ii][strainIndex_] = updating.getHapIndex(updateIndex);
      currentLLks_[ii] = updating.getNewLLKIndex(updateIndex);
      updateIndex++;

    }

    for (size_t siteI = 0; siteI < length; siteI++) {

      mcmcSample_->siteOfOneSwitchOne[start + siteI] += updating.getOneSwitchOneIndex(siteI);
      mcmcSample_->siteOfOneMissCopyOne[start + siteI] += updating.getOneMissCopyOneIndex(siteI);
      /// (K) These modified from site... to currentsite...
      mcmcSample_->currentsiteOfOneSwitchOne[start + siteI] = updating.getOneSwitchOneIndex(siteI);
      mcmcSample_->currentsiteOfOneMissCopyOne[start + siteI] = updating.getOneMissCopyOneIndex(siteI);

    }

  }

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::McmcMachinery::updatePairHaps() {

  if (kStrain() == 1) {

    return;

  }

  kgl::ExecEnv::log().vinfo(" Update Pair Hap ");

  findUpdatingStrainPair();

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts()[chromi];
    size_t length = dEploidIO_->getIndexPosition(chromi).size();

    kgl::ExecEnv::log().vinfo("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdatePairHap updating(dEploidIO_->getRefCount(),
                           dEploidIO_->getAltCount(),
                           dEploidIO_->getPlaf(),
                           currentExpectedWsaf_,
                           currentProp_,
                           currentHap_,
                           hapRg_,
                           start,
                           length,
                           panel_,
                           dEploidIO_->getMissCopyProb(),
                           dEploidIO_->scalingFactor(),
                           dEploidIO_->forbidCopyFromSame(),
                           strainIndex1_,
                           strainIndex2_);

    updating.core(dEploidIO_->getRefCount(),
                  dEploidIO_->getAltCount(),
                  dEploidIO_->getPlaf(),
                  currentExpectedWsaf_,
                  currentProp_,
                  currentHap_);

    size_t updateIndex = 0;
    for (size_t ii = start; ii < (start + length); ii++) {

      currentHap_[ii][strainIndex1_] = updating.getHap1Index(updateIndex);
      currentHap_[ii][strainIndex2_] = updating.getHap2Index(updateIndex);
      currentLLks_[ii] = updating.getNewLLKIndex(updateIndex);
      updateIndex++;

    }

    for (size_t siteI = 0; siteI < length; siteI++) {

      mcmcSample_->siteOfTwoSwitchOne[start + siteI] += updating.getTwoSwitchOneIndex(siteI);
      mcmcSample_->siteOfTwoMissCopyOne[start + siteI] += updating.getTwoMissCopyOneIndex(siteI);
      mcmcSample_->siteOfTwoSwitchTwo[start + siteI] += updating.getTwoSwitchTwoIndex(siteI);
      mcmcSample_->siteOfTwoMissCopyTwo[start + siteI] += updating.getTwoMissCopyTwoIndex(siteI);
      mcmcSample_->currentsiteOfTwoSwitchOne[start + siteI] = updating.getTwoSwitchOneIndex(siteI);
      mcmcSample_->currentsiteOfTwoMissCopyOne[start + siteI] = updating.getTwoMissCopyOneIndex(siteI);
      mcmcSample_->currentsiteOfTwoSwitchTwo[start + siteI] = updating.getTwoSwitchTwoIndex(siteI);
      mcmcSample_->currentsiteOfTwoMissCopyTwo[start + siteI] = updating.getTwoMissCopyTwoIndex(siteI);

    }

  }

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::McmcMachinery::findUpdatingStrainSingle() {

  std::vector<double> eventProb(kStrain_, 1);

  Utility::normalizeBySum(eventProb);

  strainIndex_ = Utility::sampleIndexGivenProp(mcmcEventRg_, eventProb);

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

  kgl::ExecEnv::log().vinfo("  Updating haplotype (1): {} and haplotype (2): {}", strainIndex1_, strainIndex2_);

}


