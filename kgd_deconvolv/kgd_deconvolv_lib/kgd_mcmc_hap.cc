//
// Created by kellerberrin on 25/06/18.
//


#include <random>
#include <cstdio>
#include <limits>       // std::numeric_limits< double >::min()
#include <cmath>       // ceil
#include "kgd_deconvolv_app.h"
#include "kgd_global.h"     // dout
#include "kgd_update_haplotype.h"
#include "kgd_update_single_haplotype.h"  // chromPainting
#include "kgd_update_pair_haplotype.h"
#include "kgd_mcmc_hap.h"
#include "kgd_utility.h"


namespace kgd = kellerberrin::deconvolv;


// initialiseMCMCmachinery
kgd::MCMCHAP::MCMCHAP(std::shared_ptr<DEploidIO> dEploidIO,
                                  std::shared_ptr<McmcSample> mcmcSample,
                                  std::shared_ptr<RandomGenerator> randomGenerator) : MCMCBASE(dEploidIO, mcmcSample, randomGenerator) {

  mcmcEventRg_ = randomGenerator_;
  panel_ = dEploidIO_->getPanel();

  calcMaxIteration(dEploidIO_->getMcmcSample(), dEploidIO_->getMcmcMachineryRate(), dEploidIO_->getMcmcBurn());

  MN_LOG_TITRE = 0.0;
  SD_LOG_TITRE = dEploidIO_->parameterSigma();
  PROP_SCALE = 40.0;

  initializeMcmcChain();

}


void kgd::MCMCHAP::initializeMcmcChain() {
  // Initialization

  ExecEnv::log().info("############ initializeMcmcChain() ###########");

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

  if (dEploidIO_->doAllowInbreeding()) {

    initializeUpdateReferencePanel(panel_->truePanelSize() + kStrain() - 1);

  }

  mcmcSample_->setVectorSize(nLoci());

  assert (doutProp());
  assert (doutLLK());

}




int kgd::MCMCHAP::sampleMcmcEvent() {

  int eventInt = mcmcEventRg_->sampleInt(3);

  if ((eventInt == 0) && (dEploidIO_->doUpdateProp())) {

    updateProportion();

  } else if ((eventInt == 1) && (dEploidIO_->doUpdateSingle())) {

    updateSingleHap();

  } else if ((eventInt == 2) && (dEploidIO_->doUpdatePair())) {

    updatePairHaps();

  }

  assert(doutLLK());

  return eventInt;

}


void kgd::MCMCHAP::finalizeMcmc() {


  mcmcSample_->setHap(currentHap_);

  writeLastFwdProb(false);

  dEploidIO_->setFinalProp(mcmcSample_->getProportion().back());

  mcmcSample_->divideSiteVectors(static_cast<double>(total_MCMC_iterations()));

}


void kgd::MCMCHAP::initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo) {

  if (dEploidIO_->doAllowInbreeding()) {

    return;

  }

  panel_->initializeUpdatePanel(inbreedingPanelSizeSetTo);

}


void kgd::MCMCHAP::updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain) {

  if (burnin_MCMC_iterations() > current_MCMC_iteration()) {

    return;

  }

  panel_->updatePanelWithHaps(inbreedingPanelSizeSetTo, excludedStrain, currentHap_);

}


double kgd::MCMCHAP::calcLogPriorTitre(std::vector<double> &tmpTitre) {

  //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));

  std::vector<double> tmp;

  for (auto const &value: tmpTitre) {

    tmp.push_back(std::log(Utility::normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)));

  }

  return Utility::sumOfVec(tmp);

}


void kgd::MCMCHAP::updateProportion() {

  ExecEnv::log().vinfo("MCMC Attempt of update proportion");

  if (kStrain() < 2) {

    ExecEnv::log().error("MCMC Update failed (only 1 strain)");
    return;

  }

  // calculate dt
  std::vector<double> tmpTitre = calcTmpTitre();
  std::vector<double> tmpProp = titre2prop(tmpTitre);

  if (Utility::min_value(tmpProp) < 0 || Utility::max_value(tmpProp) > 1) {

    ExecEnv::log().vinfo("MCMC Update failed, tmpProp < 0 or tmpProp > 1");
    return;

  }

  std::vector<double> tmpExpectedWsaf = calcExpectedWsaf(tmpProp);

  std::vector<double> tmpLLKs = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                                  dEploidIO_->getAltCount(),
                                                  tmpExpectedWsaf,
                                                  0,
                                                  tmpExpectedWsaf.size(),
                                                  dEploidIO_->scalingFactor());

  double diffLLKs = deltaLLKs(tmpLLKs);

  double ratio_likelihood = std::exp(diffLLKs);

  double tmpLogPriorTitre = calcLogPriorTitre(tmpTitre);

  double prior_prop_ratio = std::exp(tmpLogPriorTitre - currentLogPriorTitre_);

  double proposal_ratio = prior_prop_ratio * ratio_likelihood;

  double acceptance_draw = propRg_->sample();

  if (acceptance_draw > proposal_ratio) {

    return;  // proposal rejected.

  }


  incrementAccept();
  currentExpectedWsaf_ = tmpExpectedWsaf;
  currentLLks_ = tmpLLKs;
  currentLogPriorTitre_ = tmpLogPriorTitre;
  currentTitre_ = tmpTitre;
  currentProp_ = tmpProp;

  assert (doutProp());

}


double kgd::MCMCHAP::deltaLLKs(std::vector<double> &newLLKs) {

  std::vector<double> tmpdiff = Utility::vecDiff(newLLKs, currentLLks_);

  return Utility::sumOfVec(tmpdiff);

}


std::vector<double> kgd::MCMCHAP::calcTmpTitre() {

  std::vector<double> tmpTitre;

  for (size_t k = 0; k < kStrain(); k++) {

    double dt = deltaXnormalVariable();
    tmpTitre.push_back(currentTitre_[k] + dt);

  }

  return tmpTitre;

}


void kgd::MCMCHAP::updateSingleHap() {

  ExecEnv::log().vinfo("McmcMachinery::updateSingleHap()");

  findUpdatingStrainSingle();

  if (dEploidIO_->doAllowInbreeding()) {

    updateReferencePanel(panel_->truePanelSize() + kStrain() - 1, strainIndex_);

  }

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts()[chromi];
    size_t length = dEploidIO_->getIndexPosition(chromi).size();
    size_t kStrain = dEploidIO_->kStrain();


    ExecEnv::log().vinfo("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdateSingleHap updating(hapRg_,
                             start,
                             length,
                             kStrain,
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

      mcmcSample_->sumSiteOfOneSwitchOne(start + siteI, updating.getOneSwitchOneIndex(siteI));
      mcmcSample_->sumSiteOfOneMissCopyOne(start + siteI, updating.getOneMissCopyOneIndex(siteI));
      /// (Keller) These modified from site... to currentsite...
      mcmcSample_->setCurrentsiteOfOneSwitchOne(start + siteI, updating.getOneSwitchOneIndex(siteI));
      mcmcSample_->setCurrentsiteOfOneMissCopyOne(start + siteI, updating.getOneMissCopyOneIndex(siteI));

    }

  }

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}


void kgd::MCMCHAP::updatePairHaps() {

  if (kStrain() == 1) {

    return;

  }

  ExecEnv::log().vinfo(" Update Pair Hap ");

  findUpdatingStrainPair();

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts()[chromi];
    size_t length = dEploidIO_->getIndexPosition(chromi).size();
    size_t kStrain = dEploidIO_->kStrain();


    ExecEnv::log().vinfo("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdatePairHap updating(hapRg_,
                           start,
                           length,
                           kStrain,
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

      mcmcSample_->sumSiteOfTwoSwitchOne(start + siteI, updating.getTwoSwitchOneIndex(siteI));
      mcmcSample_->sumSiteOfTwoMissCopyOne(start + siteI, updating.getTwoMissCopyOneIndex(siteI));
      mcmcSample_->sumSiteOfTwoSwitchTwo(start + siteI, updating.getTwoSwitchTwoIndex(siteI));
      mcmcSample_->sumSiteOfTwoMissCopyTwo(start + siteI, updating.getTwoMissCopyTwoIndex(siteI));

      mcmcSample_->setCurrentsiteOfTwoSwitchOne(start + siteI, updating.getTwoSwitchOneIndex(siteI));
      mcmcSample_->setCurrentsiteOfTwoMissCopyOne(start + siteI, updating.getTwoMissCopyOneIndex(siteI));
      mcmcSample_->setCurrentsiteOfTwoSwitchTwo(start + siteI, updating.getTwoSwitchTwoIndex(siteI));
      mcmcSample_->setCurrentsiteOfTwoMissCopyTwo(start + siteI, updating.getTwoMissCopyTwoIndex(siteI));

    }

  }

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

}

void kgd::MCMCHAP::findUpdatingStrainSingle() {

  std::vector<double> eventProb(kStrain(), 1);

  Utility::normalizeBySum(eventProb);

  strainIndex_ = Utility::sampleIndexGivenProp(mcmcEventRg_, eventProb);

  ExecEnv::log().vinfo("Updating haplotype (1 of 1): {}", strainIndex_);

}


void kgd::MCMCHAP::findUpdatingStrainPair() {

  std::vector<size_t> strainIndex(2, 0);

  int t = 0; // total input records dealt with
  int m = 0; // number of items selected so far
  double u;

  while (m < 2) {

    u = mcmcEventRg_->sample(); // call a uniform(0,1) kgd_random number generator

    if ((kStrain() - t) * u < 2 - m) {

      strainIndex[m] = t;
      m++;

    }

    t++;

  }

  strainIndex1_ = strainIndex[0];
  strainIndex2_ = strainIndex[1];

  assert(strainIndex1_ != strainIndex2_);

  ExecEnv::log().vinfo("Updating haplotype (1 of 2): {} and haplotype (2 of 2): {}", strainIndex1_, strainIndex2_);

}


void kgd::MCMCHAP::writeLastFwdProb(bool useIBD) {

  if (not dEploidIO_->doExportPostProb()) {

    return;

  }

  for (size_t tmpk = 0; tmpk < kStrain(); tmpk++) {

    if (dEploidIO_->doAllowInbreeding()) {

      updateReferencePanel(panel_->truePanelSize() + kStrain() - 1, tmpk);

    }

    for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

      size_t start = dEploidIO_->indexOfChromStarts()[chromi];
      size_t length = dEploidIO_->getIndexPosition(chromi).size();
      size_t kStrain = dEploidIO_->kStrain();

      UpdateSingleHap updatingSingle(hapRg_,
                                     start,
                                     length,
                                     kStrain,
                                     panel_,
                                     dEploidIO_->getMissCopyProb(),
                                     dEploidIO_->scalingFactor(),
                                     tmpk);

      if (dEploidIO_->doAllowInbreeding()) {

        updatingSingle.setPanelSize(panel_->inbreedingPanelSize());

      }

      updatingSingle.core(dEploidIO_->getRefCount(),
                          dEploidIO_->getAltCount(),
                          dEploidIO_->getPlaf(),
                          currentExpectedWsaf_,
                          currentProp_,
                          currentHap_);

      dEploidIO_->writeLastSingleFwdProb(updatingSingle.getFwdProbs(), chromi, tmpk, useIBD);

    }

  }

}

