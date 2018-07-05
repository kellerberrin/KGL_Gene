//
// Created by kellerberrin on 25/06/18.
//


#include <random>
#include <cstdio>
#include <limits>       // std::numeric_limits< double >::min()
#include <cmath>       // ceil
#include "kgd_deconvolv_app.h"
#include "kgd_global.h"     // dout
#include "kgd_mcmc_ibd.h"
#include "kgd_utility.h"

namespace kgd = kellerberrin::deconvolv;


// initialiseMCMCmachinery
kgd::MCMCIBD::MCMCIBD(std::shared_ptr<DEploidIO> dEploidIO,
                      std::shared_ptr<McmcSample> mcmcSample,
                      std::shared_ptr<RandomGenerator> randomGenerator)
  : MCMCBASE(dEploidIO, mcmcSample, randomGenerator, dEploidIO->ibdParameters()) {

  initializeMcmcChain();

}




void kgd::MCMCIBD::initializeMcmcChain() {
  // Initialization

  ExecEnv::log().info("############ initializeMcmcChain() ###########");

  initializeHap();
  initializeExpectedWsaf(titre_proportions_.Proportions()); // This requires currentHap_ and currentProp_
  currentLLks_ = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                   dEploidIO_->getAltCount(),
                                   currentExpectedWsaf_,
                                   0,
                                   currentExpectedWsaf_.size(),
                                   dEploidIO_->ibdParameters().proposalUpdateScaling());

  ibdInitializeEssentials(MCMCParameters_.baseCountError());

  mcmcSample_->setVectorSize(nLoci());

}


int kgd::MCMCIBD::sampleMcmcEvent() {

  ibdSampleMcmcEventStep();

  return 0;

}


void kgd::MCMCIBD::finalizeMcmc() {

  mcmcSample_->setHap(currentHap_);

  dEploidIO_->setFinalProp(mcmcSample_->getProportion().back());

  mcmcSample_->divideSiteVectors(static_cast<double>(total_MCMC_iterations()));

  for (size_t atSiteI = 0; atSiteI < nLoci(); atSiteI++) {

    ibdPath.IBDPathChangeAt(atSiteI,  static_cast<double>(total_MCMC_iterations()));

  }

  ExecEnv::log().info("Proportion update acceptance rate: {}", acceptCount() / (kStrain() * 1.0 * total_MCMC_iterations()));

  dEploidIO_->setInitialProp(averageProportion(mcmcSample_->getProportion()));
  dEploidIO_->setInitialPropWasGiven(true);
  dEploidIO_->setDoUpdateProp(false);
  dEploidIO_->setInitialHap(mcmcSample_->getHap());
  dEploidIO_->setInitialHapWasGiven(true);

}


std::vector<double> kgd::MCMCIBD::averageProportion(const std::vector<std::vector<double> > &proportion) {

  assert(proportion.size() > 0);

  std::vector<double> average_proportion;

  for (size_t i = 0; i < kStrain(); i++) {

    double strain_average = 0;

    for (auto prop_strain : proportion) {

      strain_average += prop_strain[i];

    }

    strain_average = strain_average / static_cast<double>(proportion.size());

    average_proportion.push_back(strain_average);

  }

  Utility::normalizeBySum(average_proportion);

  return average_proportion;

}


void kgd::MCMCIBD::ibdInitializeEssentials(double err) {

  ibdPath.init(*dEploidIO_, hapRg_);

  std::vector<double> llkOfData;

  for (size_t i = 0; i < nLoci(); i++) {

    double wsaf = dEploidIO_->getAltCount()[i] / (dEploidIO_->getRefCount()[i] + dEploidIO_->getAltCount()[i] + 0.00000000000001);

    double adjustedWsaf = (wsaf * (1 - err)) + ((1 - wsaf) * err);

    double loci_likelihood = Utility::logBetaPdf(adjustedWsaf, ibdPath.getLogLikelihoodSurface()[i][0], ibdPath.getLogLikelihoodSurface()[i][1]);

    llkOfData.push_back(loci_likelihood);

  }

  ExecEnv::log().info("LLK of data = {}", Utility::sumOfVec(llkOfData));

}


void kgd::MCMCIBD::ibdSampleMcmcEventStep() {

  // Update the idb path.
  ibdPath.McmcUpdateStep(titre_proportions_.Proportions());

  //#Get haplotypes and update LLK for each site
  ibdUpdateHaplotypesFromPrior();

  ///#Given current haplotypes, sample titres 1 by 1 using MH
  ibdUpdateProportionGivenHap();

  currentExpectedWsaf_ = calcExpectedWsaf(titre_proportions_.Proportions());

}


void kgd::MCMCIBD::ibdUpdateHaplotypesFromPrior() {

  for (size_t loci = 0; loci < nLoci(); ++loci) {

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      currentHap_[loci][strain] = ibdPath.UpdateHaplotypesFromPrior(strain, loci);

    }

  }

}


void kgd::MCMCIBD::ibdUpdateProportionGivenHap() {

  // For each strain perform a Metropolis-Hastings update of currentTitre, currentProportion
  std::vector<double> llkAtAllSites = computeLlkAtAllSites(titre_proportions_.Proportions(), MCMCParameters_.baseCountError());

  for (size_t i = 0; i < kStrain(); i++) {

    MCMCTITRE proposal(titre_proportions_);

    proposal.updateTitreIndex(i);

    std::vector<double> llk_loci_vector = computeLlkAtAllSites(proposal.Proportions(), MCMCParameters_.baseCountError());

    double prior_titre_ratio = proposal.calcPriorTitreIndex(i) / titre_proportions_.calcPriorTitreIndex(i);
    double likelihood_ratio = exp(Utility::sumOfVec(llk_loci_vector) - Utility::sumOfVec(llkAtAllSites));
    double proposal_ratio = prior_titre_ratio * likelihood_ratio;

    double acceptance_draw = propRg_->sample();

    ExecEnv::log().info("{}, IDB [{}]: {}/{}, h_ratio: {}, l_ratio: {}, h*r: {}, accept: {}, {}",
                        current_MCMC_iteration(), i, proposal.Proportions()[i], titre_proportions_.Proportions()[i],
                        prior_titre_ratio, likelihood_ratio, proposal_ratio, acceptance_draw, acceptance_draw <= proposal_ratio ? "ACCEPT" : "REJECT");

    if (acceptance_draw <= proposal_ratio) {
      // Accept
      titre_proportions_ = proposal;
      llkAtAllSites = llk_loci_vector;
      incrementAccept();

    }

  }

  currentLLks_ = llkAtAllSites;

}


std::vector<double> kgd::MCMCIBD::computeLlkAtAllSites(const std::vector<double>& proportion, double read_error_prob) {

  std::vector<double> llk_vector;

  for (size_t loci = 0; loci < nLoci(); ++loci) {

    double loci_proportion = 0;

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      loci_proportion += (double) currentHap_[loci][strain] * proportion[strain];

    }

    double loci_prop_err = (loci_proportion * (1 - read_error_prob)) + ((1 - loci_proportion) * read_error_prob);

    double llk_loci = Utility::logBetaPdf(loci_prop_err, ibdPath.getLogLikelihoodSurface()[loci][0], ibdPath.getLogLikelihoodSurface()[loci][1]);

    llk_vector.push_back(llk_loci);

  }

  return llk_vector;

}


void kgd::MCMCIBD::recordMcmcMachinery() {

  mcmcSample_->addProportion(titre_proportions_.Proportions());
  mcmcSample_->addSumLLKs(Utility::sumOfVec(currentLLks_));
  mcmcSample_->addMove(eventType());

  // Accumulate expectedWSAF for computing the mean expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); ++i) {

    cumExpectedWsaf_[i] += currentExpectedWsaf_[i];

  }

}

