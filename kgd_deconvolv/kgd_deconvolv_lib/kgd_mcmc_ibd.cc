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

//#define MCMCIBD_DEBUG 1

// initialiseMCMCmachinery
kgd::MCMCIBD::MCMCIBD(std::shared_ptr<DEploidIO> dEploidIO,
                      std::shared_ptr<McmcSample> mcmcSample)
  : MCMCBASE(dEploidIO, mcmcSample, dEploidIO->ibdParameters()) {

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

  theta_accept_ = 0;

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

    ibdPath_.IBDPathChangeAt(atSiteI,  static_cast<double>(total_MCMC_iterations()));

  }

//  ExecEnv::log().info("Proportion update acceptance rate: {}", acceptCount() / (kStrain() * 1.0 * total_MCMC_iterations()));
  ExecEnv::log().info("Proportion update acceptance rate: {}", acceptCount() / (1.0 * total_MCMC_iterations()));

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

  ibdPath_.init(*dEploidIO_);

  std::vector<double> llkOfData;

  for (size_t site = 0; site < nLoci(); ++site) {

    double wsaf = dEploidIO_->getAltCount()[site] / (dEploidIO_->getRefCount()[site] + dEploidIO_->getAltCount()[site] + 0.00000000000001);

    double adjustedWsaf = (wsaf * (1 - err)) + ((1 - wsaf) * err);

    double loci_likelihood = ibdPath_.siteLogBetaLLK(site, adjustedWsaf);

    llkOfData.push_back(loci_likelihood);

  }

  ExecEnv::log().info("LLK of data = {}", Utility::sumOfVec(llkOfData));

}


void kgd::MCMCIBD::ibdSampleMcmcEventStep() {

  static BetaDistribution random_beta(2.0,2.0);
  //  const double theta_proposal_size = 0.1;
  std::vector<double> llkAtAllSites = computeLlkAtAllSites(titre_proportions_.Proportions(), MCMCParameters_.baseCountError());

  double previous_theta = ibdPath_.theta();

//  double proposed_theta = previous_theta;

//  do {

//    proposed_theta += ((propRg_->sample() * theta_proposal_size) - (theta_proposal_size / 2.0));

//  } while (proposed_theta < 0 or proposed_theta > 1);

//  double theta = random_unit_.generate(entropy_source_.generator());
  double theta = random_beta.random(entropy_source_.generator());
  ibdPath_.setTheta(theta);

  // Update the idb path.
  ibdPath_.McmcUpdateStep(titre_proportions_.Proportions());

  //#Get haplotypes and update LLK for each site
  ibdUpdateHaplotypesFromPrior();

  /// Experimental: update theta using MCMC
  ibdUpdateThetaGivenHap(previous_theta, llkAtAllSites);

  ///#Given current haplotypes, sample titres 1 by 1 using MH
  updateProportion();
//  ibdUpdateProportionGivenHap();

  currentExpectedWsaf_ = calcExpectedWsaf(titre_proportions_.Proportions());

}


void kgd::MCMCIBD::ibdUpdateHaplotypesFromPrior() {

  for (size_t site_i = 0; site_i < nLoci(); ++site_i) {

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      currentHap_[site_i][strain] = ibdPath_.UpdateHaplotypesFromPrior(strain, site_i);

    }

  }

}


void kgd::MCMCIBD::ibdUpdateThetaGivenHap(double previous_theta, const std::vector<double>& llkAtAllSites) {


  std::vector<double> updated_llkAtAllSites = computeLlkAtAllSites(titre_proportions_.Proportions(), MCMCParameters_.baseCountError());

  double prior_ratio = 1.0;  // Since we are sampling from the uniform distribution.

  double log_likelihood_ratio = Utility::sumOfVec(updated_llkAtAllSites) - Utility::sumOfVec(llkAtAllSites);
  double likelihood_ratio = std::exp(log_likelihood_ratio);
  double proposal_ratio = prior_ratio * likelihood_ratio;
  double acceptance_draw = random_unit_.random(entropy_source_.generator());

#ifdef MCMCIBD_DEBUG

  ExecEnv::log().info("{}, Theta: {}/{}, ll_r: {}, l_r: {}, h*r: {}, accept: {}, ratio: {}, {}",
                      current_MCMC_iteration(), ibdPath_.theta(), previous_theta, log_likelihood_ratio,
                      likelihood_ratio, proposal_ratio, acceptance_draw,
                      static_cast<double>(theta_accept_)/static_cast<double>(current_MCMC_iteration()),
                      acceptance_draw <= proposal_ratio ? "ACCEPT" : "REJECT");

#endif

  if (acceptance_draw > proposal_ratio) {
    // Reject
    ibdPath_.setTheta(previous_theta);

  } else {

    theta_accept_ += 1;

  }

}



void kgd::MCMCIBD::updateProportion() {


  // For each strain perform a Metropolis-Hastings update of currentTitre, currentProportion
  std::vector<double> llkAtAllSites = computeLlkAtAllSites(titre_proportions_.Proportions(), MCMCParameters_.baseCountError());

  MCMCTITRE proposal(titre_proportions_);

  proposal.updateTitre();

  std::vector<double> llk_loci_vector = computeLlkAtAllSites(proposal.Proportions(), MCMCParameters_.baseCountError());

  double prior_prop_ratio = std::exp(proposal.calcLogPriorTitre() - titre_proportions_.calcLogPriorTitre());

  double log_likelihood_ratio = Utility::sumOfVec(llk_loci_vector) - Utility::sumOfVec(llkAtAllSites);

  double likelihood_ratio = std::exp(log_likelihood_ratio);

  double proposal_ratio = prior_prop_ratio * likelihood_ratio * proposal.hastingsRatio();

  double acceptance_draw = random_unit_.random(entropy_source_.generator());

#ifdef MCMCIBD_DEBUG

  double ratio = static_cast<double>(acceptCount()) / static_cast<double>(current_MCMC_iteration());

  ExecEnv::log().info("{}, Update: {}, p_r: {}, ll_r: {}, l_r: {}, p*l: {}, accept: {}, {}, ratio: {}",
                      current_MCMC_iteration(), proposal.proportionsText(),
                      prior_prop_ratio, log_likelihood_ratio, likelihood_ratio, proposal_ratio, acceptance_draw,
                      acceptance_draw <= proposal_ratio ? "ACCEPT" : "REJECT", ratio);

#endif

  if (acceptance_draw <= proposal_ratio) {
  // Accept
    incrementAccept();
    llkAtAllSites = llk_loci_vector;
    titre_proportions_ = proposal;

  }

}



void kgd::MCMCIBD::ibdUpdateProportionGivenHap() {

  // For each strain perform a Metropolis-Hastings update of currentTitre, currentProportion
  std::vector<double> llkAtAllSites = computeLlkAtAllSites(titre_proportions_.Proportions(), MCMCParameters_.baseCountError());

  for (size_t k = 0; k < kStrain(); k++) {

    MCMCTITRE proposal(titre_proportions_);

    proposal.updateTitreIndex(k);

    std::vector<double> llk_loci_vector = computeLlkAtAllSites(proposal.Proportions(), MCMCParameters_.baseCountError());

    double prior_titre_ratio = proposal.calcPriorTitreIndex(k) / titre_proportions_.calcPriorTitreIndex(k);
    double log_likelihood_ratio = Utility::sumOfVec(llk_loci_vector) - Utility::sumOfVec(llkAtAllSites);
    double likelihood_ratio = std::exp(log_likelihood_ratio);
    double proposal_ratio = prior_titre_ratio * likelihood_ratio;
    double acceptance_draw = random_unit_.random(entropy_source_.generator());

#ifdef MCMCIBD_DEBUG

    ExecEnv::log().info("{}, IDB [{}]: {}/{}, h_ratio: {}, ll_r: {}, l_r: {}, h*r: {}, accept: {}, {}",
                        current_MCMC_iteration(), k, proposal.Proportions()[k], titre_proportions_.Proportions()[k],
                        prior_titre_ratio, log_likelihood_ratio, likelihood_ratio, proposal_ratio, acceptance_draw,
                        acceptance_draw <= proposal_ratio ? "ACCEPT" : "REJECT");

#endif

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

      loci_proportion += currentHap_[loci][strain] * proportion[strain];

    }

    double loci_prop_err = (loci_proportion * (1 - read_error_prob)) + ((1 - loci_proportion) * read_error_prob);

    double llk_loci = ibdPath_.siteLogBetaLLK(loci, loci_prop_err);

    llk_vector.push_back(llk_loci);

  }

  return llk_vector;

}


void kgd::MCMCIBD::recordMcmcMachinery() {

  mcmcSample_->addProportion(titre_proportions_.Proportions());
  double sumLLKS = Utility::sumOfVec(currentLLks_);
  mcmcSample_->addSumLLKs(sumLLKS);
  mcmcSample_->addMove(eventType());

  ExecEnv::log().info("Iterations: {}, Proportions: {} Sum Log Likelihoods: {}",
                           current_MCMC_iteration(), titre_proportions_.proportionsText(), sumLLKS);

  // Accumulate expectedWSAF for computing the mean expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); ++i) {

    cumExpectedWsaf_[i] += currentExpectedWsaf_[i];

  }

}

