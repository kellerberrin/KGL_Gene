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
  : MCMCBASE(dEploidIO, mcmcSample, randomGenerator),
    titre_proportions_(kStrain(), 0.0 /* mean */, dEploidIO_->ibdSigma(), 40.0 /*update */) {

  calcMaxIteration(100, 10, 0.5);  // TODO: Get these values from Deploid_io.

  MN_LOG_TITRE = 0.0;
  SD_LOG_TITRE = dEploidIO_->ibdSigma();
  PROP_SCALE = 40.0;

  initializeMcmcChain();

}




void kgd::MCMCIBD::initializeMcmcChain() {
  // Initialization

  ExecEnv::log().info("############ initializeMcmcChain() ###########");

  initializeTitre();
  initializeHap();
  initializeProp();


  initializeExpectedWsaf(); // This requires currentHap_ and currentProp_
  currentLLks_ = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                   dEploidIO_->getAltCount(),
                                   currentExpectedWsaf_,
                                   0,
                                   currentExpectedWsaf_.size(),
                                   dEploidIO_->scalingFactor());

  ibdInitializeEssentials();

  mcmcSample_->setVectorSize(nLoci());

  assert (doutProp());
  assert (doutLLK());

}


int kgd::MCMCIBD::sampleMcmcEvent() {

  ibdSampleMcmcEventStep();

  assert(doutProp());
  assert(doutLLK());

  return 0;

}


void kgd::MCMCIBD::finalizeMcmc() {

  mcmcSample_->setHap(currentHap_);

//  writeLastFwdProb();

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


void kgd::MCMCIBD::initializePropIBD() {


  if (dEploidIO_->initialPropWasGiven()) {

    titre_proportions_.proportion2Titre(dEploidIO_->getInitialProp());

  }

  currentProp_ = (dEploidIO_->initialPropWasGiven()) ? dEploidIO_->getInitialProp() : titre2prop(currentTitre_);

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

  initializePropIBD();

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
  ibdPath.McmcUpdateStep(currentProp_);

  //#Get haplotypes and update LLK for each site
  ibdUpdateHaplotypesFromPrior();

  ////#Given current haplotypes, sample titres 1 by 1 using MH

  ibdUpdateProportionGivenHap();

  currentExpectedWsaf_ = calcExpectedWsaf(currentProp_);

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

  std::vector<double> llkAtAllSites = computeLlkAtAllSites(currentProp_);

  for (size_t i = 0; i < kStrain(); i++) {

    MCMCTITRE proposal(titre_proportions_);

    proposal.updateTitreIndex(i);

    double titre_i = currentTitre_[i];

    std::vector<double> oldProp = currentProp_;

    currentTitre_[i] += ((stdNorm_->genReal() * (SD_LOG_TITRE / PROP_SCALE)) + MN_LOG_TITRE); // tit.0[i]+rnorm(1, 0, scale.t.prop);
    currentProp_ = titre2prop(currentTitre_);

    std::vector<double> llk_loci_vector = computeLlkAtAllSites(currentProp_);

//    double hastings_ratio = Utility::normal_pdf(currentTitre_[i], 0, 1) / Utility::normal_pdf(titre_i, 0, 1);
    double prior_ratio = Utility::normal_pdf(currentTitre_[i], MN_LOG_TITRE, SD_LOG_TITRE) / Utility::normal_pdf(titre_i, MN_LOG_TITRE, SD_LOG_TITRE);
    double likelihood_ratio = exp(Utility::sumOfVec(llk_loci_vector) - Utility::sumOfVec(llkAtAllSites));
    double proposal_ratio = prior_ratio * likelihood_ratio;

    double acceptance_draw = propRg_->sample();

    ExecEnv::log().info("IDB prop[{}]: {}, h_ratio: {}, l_ratio: {}, pr_ratio: {}, accept_draw: {}, {}",
                        i, currentProp_[i], prior_ratio, likelihood_ratio, proposal_ratio, acceptance_draw, acceptance_draw <= proposal_ratio ? "ACCEPT" : "REJECT");

    if (acceptance_draw <= proposal_ratio) {
      // Accept
      llkAtAllSites = llk_loci_vector;
      incrementAccept();

    } else {
      // Reject
      currentTitre_[i] = titre_i;
      currentProp_ = oldProp;

    }

  }

  currentLLks_ = llkAtAllSites;

}


std::vector<double> kgd::MCMCIBD::computeLlkAtAllSites(const std::vector<double>& proportion, double err) {

  std::vector<double> llk_vector;

  for (size_t loci = 0; loci < nLoci(); ++loci) {

    double loci_proportion = 0;

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      loci_proportion += (double) currentHap_[loci][strain] * proportion[strain];

    }

    double loc_proportion_err = (loci_proportion * (1 - err)) + ((1 - loci_proportion) * err);

    double llk_loci = Utility::logBetaPdf(loc_proportion_err, ibdPath.getLogLikelihoodSurface()[loci][0], ibdPath.getLogLikelihoodSurface()[loci][1]);

    llk_vector.push_back(llk_loci);

  }

  return llk_vector;

}


