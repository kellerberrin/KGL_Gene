//
// Created by kellerberrin on 25/06/18.
//


#include <random>
#include <stdio.h>
#include <limits>       // std::numeric_limits< double >::min()
#include <math.h>       // ceil
#include "kgd_deconvolv_app.h"
#include "kgd_global.h"     // dout
#include "kgd_mcmc_ibd.h"
#include "kgd_utility.h"

namespace kgd = kellerberrin::deconvolv;


// initialiseMCMCmachinery
kgd::MCMCIBD::MCMCIBD(std::shared_ptr<DEploidIO> dEploidIO,
                      std::shared_ptr<McmcSample> mcmcSample,
                      std::shared_ptr<RandomGenerator> randomGenerator) : MCMCBASE(dEploidIO, mcmcSample, randomGenerator) {

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

  if (dEploidIO_->doAllowInbreeding()) {

    initializeUpdateReferencePanel(panel_->truePanelSize() + kStrain() - 1);

  }

  ibdInitializeEssentials();

  mcmcSample_->setVectorSize(nLoci());

  assert (doutProp());
  assert (doutLLK());

}



void kgd::MCMCIBD::initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo) {

  if (dEploidIO_->doAllowInbreeding()) {

    return;

  }

  panel_->initializeUpdatePanel(inbreedingPanelSizeSetTo);

}



double kgd::MCMCIBD::calcLogPriorTitre(std::vector<double> &tmpTitre) {
  //sum(dnorm(titre, MN_LOG_TITRE, SD_LOG_TITRE, log=TRUE));

  std::vector<double> tmp;

  for (auto const &value: tmpTitre) {

    tmp.push_back(log(Utility::normal_pdf(value, MN_LOG_TITRE, SD_LOG_TITRE)));

  }

  return Utility::sumOfVec(tmp);

}



void kgd::MCMCIBD::initializePropIBD() {

  currentProp_ = (dEploidIO_->initialPropWasGiven()) ? dEploidIO_->getInitialProp() : titre2prop(currentTitre_);

}


void kgd::MCMCIBD::finalizeMcmc() {

  mcmcSample_->setHap(currentHap_);

//  writeLastFwdProb();

  dEploidIO_->setFinalProp(mcmcSample_->getProportion().back());

  mcmcSample_->divideSiteVectors(static_cast<double>(maxIteration_));

  dEploidIO_->writeMcmcRelated(mcmcSample_, true);

  for (size_t atSiteI = 0; atSiteI < nLoci(); atSiteI++) {

    ibdPath.IBDPathChangeAt(atSiteI,  static_cast<double>(maxIteration_));

  }

  ExecEnv::log().info("Proportion update acceptance rate: {}", acceptCount() / (kStrain() * 1.0 * maxIteration_));

  dEploidIO_->setInitialProp(averageProportion(mcmcSample_->getProportion()));
  dEploidIO_->setInitialPropWasGiven(true);
  dEploidIO_->setDoUpdateProp(false);
  dEploidIO_->setInitialHap(mcmcSample_->getHap());
  dEploidIO_->setInitialHapWasGiven(true);

}


std::vector<double> kgd::MCMCIBD::averageProportion(const std::vector<std::vector<double> > &proportion) {

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


int kgd::MCMCIBD::sampleMcmcEvent() {

  recordingMcmcBool_ = (currentMcmcIteration_ > mcmcThresh_ && currentMcmcIteration_ % McmcMachineryRate_ == 0);

  ibdSampleMcmcEventStep();

  assert(doutProp());
  assert(doutLLK());

  if (recordingMcmcBool_) {

    recordMcmcMachinery();

  }

  return 0;

}


void kgd::MCMCIBD::ibdInitializeEssentials() {

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


void kgd::MCMCIBD::ibdSampleMcmcEventStep() {

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


void kgd::MCMCIBD::ibdUpdateHaplotypesFromPrior() {

  for (size_t loci = 0; loci < nLoci(); ++loci) {

    for (size_t strain = 0; strain < kStrain(); ++strain) {

      currentHap_[loci][strain] = ibdPath.UpdateHaplotypesFromPrior(strain, loci);

    }

  }

}


void kgd::MCMCIBD::ibdUpdateProportionGivenHap(std::vector<double> &llkAtAllSites) {

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
      incrementAccept();

    } else {

      currentTitre_[i] = v0;
      currentProp_ = oldProp;

    }

  }

}


std::vector<double> kgd::MCMCIBD::computeLlkAtAllSites(double err) {

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


