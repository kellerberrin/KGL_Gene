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


//#define MCMCHAP_DEBUG 1

namespace kgd = kellerberrin::deconvolv;


// initialiseMCMCmachinery
kgd::MCMCHAP::MCMCHAP(std::shared_ptr<DEploidIO> dEploidIO,
                      std::shared_ptr<McmcSample> mcmcSample)
  : MCMCBASE(dEploidIO, mcmcSample, dEploidIO->hapParameters()) , rand_strain_(dEploidIO_->kStrain()) {

  panel_ = dEploidIO_->getPanel();

  initializeMcmcChain();

}


void kgd::MCMCHAP::initializeMcmcChain() {
  // Initialization

  ExecEnv::log().info("############ initializeMcmcChain() ###########");

  initializeHap();
  initializeExpectedWsaf(titre_proportions_.Proportions()); // This requires currentHap_ and currentProp_
  currentLLks_ = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                   dEploidIO_->getAltCount(),
                                   currentExpectedWsaf_,
                                   0,
                                   currentExpectedWsaf_.size(),
                                   MCMCParameters_.proposalUpdateScaling());

  if (dEploidIO_->doAllowInbreeding()) {

    initializeUpdateReferencePanel(panel_->truePanelSize() + kStrain() - 1);

  }

  mcmcSample_->setVectorSize(nLoci());

}


int kgd::MCMCHAP::sampleMcmcEvent() {

  HapUpdateType Update = rand_update_type_.generateUpdate();

  switch(Update) {

    case HapUpdateType::UPDATE_PROPORTION:
      if (dEploidIO_->doUpdateProp()) updateProportion();
      break;

    case HapUpdateType::UPDATE_SINGLE_HAPLOTYPE:
      if (dEploidIO_->doUpdateSingle()) updateSingleHap(titre_proportions_.Proportions());
      break;

    case HapUpdateType::UPDATE_PAIR_HAPLOTYPE:
      if (dEploidIO_->doUpdatePair()) updatePairHaps(titre_proportions_.Proportions());
      break;

  }

  return static_cast<int>(Update);

}


void kgd::MCMCHAP::finalizeMcmc() {


  mcmcSample_->setHap(currentHap_);

  writeLastFwdProb(titre_proportions_.Proportions(), false);

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


void kgd::MCMCHAP::updateProportion() {


  MCMCTITRE proposal(titre_proportions_);

  proposal.updateTitre();

  std::vector<double> tmpExpectedWsaf = calcExpectedWsaf(proposal.Proportions());

  std::vector<double> tmpLLKs = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                                  dEploidIO_->getAltCount(),
                                                  tmpExpectedWsaf,
                                                  0,
                                                  tmpExpectedWsaf.size(),
                                                  MCMCParameters_.proposalUpdateScaling());

  double diffLLKs = deltaLLKs(tmpLLKs);

  double ratio_likelihood = std::exp(diffLLKs);

  double prior_prop_ratio = std::exp(proposal.calcLogPriorTitre() - titre_proportions_.calcLogPriorTitre());

  double proposal_ratio = prior_prop_ratio * ratio_likelihood * proposal.hastingsRatio();

  double acceptance_draw = random_unit_.random(entropy_source_.generator());

#ifdef MCMCHAP_DEBUG

  ExecEnv::log().info("{}, Update: {}, p_ratio: {}, l_ratio: {}, p*l: {}, accept: {}, {}",
                      current_MCMC_iteration(), proposal.proportionsText(),
                      prior_prop_ratio, ratio_likelihood, proposal_ratio, acceptance_draw, acceptance_draw <= proposal_ratio ? "ACCEPT" : "REJECT");

#endif

  if (acceptance_draw > proposal_ratio) {

    return;  // proposal rejected.

  }

  incrementAccept();
  currentExpectedWsaf_ = tmpExpectedWsaf;
  currentLLks_ = tmpLLKs;
  titre_proportions_ = proposal;

}


double kgd::MCMCHAP::deltaLLKs(std::vector<double> &newLLKs) {

  std::vector<double> tmpdiff = Utility::vecDiff(newLLKs, currentLLks_);

  return Utility::sumOfVec(tmpdiff);

}


void kgd::MCMCHAP::updateSingleHap(const std::vector<double>& proportions) {

  ExecEnv::log().vinfo("McmcMachinery::updateSingleHap()");

  size_t strain_index = rand_strain_.getRandomStrain();

  if (dEploidIO_->doAllowInbreeding()) {

    updateReferencePanel(panel_->truePanelSize() + kStrain() - 1, strain_index);

  }

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts()[chromi];
    size_t length = dEploidIO_->getIndexPosition(chromi).size();
    size_t kStrain = dEploidIO_->kStrain();


    ExecEnv::log().vinfo("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdateSingleHap updating(start,
                             length,
                             kStrain,
                             panel_,
                             MCMCParameters_.getMissCopyProb(),
                             MCMCParameters_.proposalUpdateScaling(),
                             strain_index);

    if (dEploidIO_->doAllowInbreeding()) {

      updating.setPanelSize(panel_->inbreedingPanelSize());

    }

    updating.core(dEploidIO_->getRefCount(),
                  dEploidIO_->getAltCount(),
                  dEploidIO_->getPlaf(),
                  currentExpectedWsaf_,
                  proportions,
                  currentHap_);

    size_t updateIndex = 0;

    for (size_t ii = start; ii < (start + length); ii++) {

      currentHap_[ii][strain_index] = updating.getHapIndex(updateIndex);
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

  currentExpectedWsaf_ = calcExpectedWsaf(proportions);

}


void kgd::MCMCHAP::updatePairHaps(const std::vector<double>& proportions) {

  if (kStrain() == 1) {

    return;

  }

  ExecEnv::log().vinfo(" Update Pair Hap ");

//  findUpdatingStrainPair();

  std::pair<size_t, size_t> strain_pair = rand_strain_.getRandomStrainPair();

  for (size_t chromi = 0; chromi < dEploidIO_->indexOfChromStarts().size(); chromi++) {

    size_t start = dEploidIO_->indexOfChromStarts()[chromi];
    size_t length = dEploidIO_->getIndexPosition(chromi).size();
    size_t kStrain = dEploidIO_->kStrain();


    ExecEnv::log().vinfo("Update Chrom with index: {}, starts at: {}, with: {} sites", chromi, start, length);

    UpdatePairHap updating(start,
                           length,
                           kStrain,
                           panel_,
                           MCMCParameters_.getMissCopyProb(),
                           MCMCParameters_.proposalUpdateScaling(),
                           dEploidIO_->forbidCopyFromSame(),
                           strain_pair.first,
                           strain_pair.second);

    updating.core(dEploidIO_->getRefCount(),
                  dEploidIO_->getAltCount(),
                  dEploidIO_->getPlaf(),
                  currentExpectedWsaf_,
                  proportions,
                  currentHap_);

    size_t updateIndex = 0;
    for (size_t ii = start; ii < (start + length); ii++) {

      currentHap_[ii][strain_pair.first] = updating.getHap1Index(updateIndex);
      currentHap_[ii][strain_pair.second] = updating.getHap2Index(updateIndex);
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

  currentExpectedWsaf_ = calcExpectedWsaf(proportions);

}


void kgd::MCMCHAP::writeLastFwdProb(const std::vector<double>& proportions, bool useIBD) {

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

      UpdateSingleHap updatingSingle(start,
                                     length,
                                     kStrain,
                                     panel_,
                                     MCMCParameters_.getMissCopyProb(),
                                     MCMCParameters_.proposalUpdateScaling(),
                                     tmpk);

      if (dEploidIO_->doAllowInbreeding()) {

        updatingSingle.setPanelSize(panel_->inbreedingPanelSize());

      }

      updatingSingle.core(dEploidIO_->getRefCount(),
                          dEploidIO_->getAltCount(),
                          dEploidIO_->getPlaf(),
                          currentExpectedWsaf_,
                          proportions,
                          currentHap_);

      dEploidIO_->writeLastSingleFwdProb(updatingSingle.getFwdProbs(), chromi, tmpk, useIBD);

    }

  }

}


void kgd::MCMCHAP::recordMcmcMachinery() {

  mcmcSample_->addProportion(titre_proportions_.Proportions());
  mcmcSample_->addSumLLKs(Utility::sumOfVec(currentLLks_));
  mcmcSample_->addMove(eventType());

  // Accumulate expectedWSAF for computing the mean expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); ++i) {

    cumExpectedWsaf_[i] += currentExpectedWsaf_[i];

  }

}


