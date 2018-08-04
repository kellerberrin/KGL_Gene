//
// Created by kellerberrin on 25/06/18.
//


#include "kgd_deconvolv_app.h"
#include "kgd_global.h"     // dout
#include "kgd_utility.h"
#include "kgd_mcmc_base.h"


namespace kgd = kellerberrin::deconvolv;


kgd::MCMCBASE::MCMCBASE(std::shared_ptr<DEploidIO> dEploidIO,
                        std::shared_ptr<McmcSample> mcmcSample,
                        const MCMCParameterObj& MCMCParameters)
  : MCMCVIRTUAL(dEploidIO, mcmcSample),
    MCMCParameters_(MCMCParameters),
    titre_proportions_(dEploidIO->kStrain(),
                       MCMCParameters.proposalMean(),
                       MCMCParameters.proposalSigma(),
                       MCMCParameters.proposalUpdateScaling()) {

  calcMaxIteration(MCMCParameters.McmcSample(),
                   MCMCParameters.McmcMachineryRate(),
                   MCMCParameters.McmcBurn());

  if (dEploidIO_->initialPropWasGiven()) {

    titre_proportions_.proportion2Titre(dEploidIO_->getInitialProp());

  }

  setKstrain(dEploidIO_->kStrain());
  setNLoci(dEploidIO_->getPlaf().size());

}


void kgd::MCMCBASE::initializeHap() {

  assert(currentHap_.size() == 0);

  if (dEploidIO_->initialHapWasGiven()) {

    currentHap_ = dEploidIO_->getInitialHap();

  } else {

    for (size_t i = 0; i < dEploidIO_->getPlaf().size(); i++) {

      double currentPlaf = dEploidIO_->getPlaf()[i];
      std::vector<double> tmpVec;

      for (size_t k = 0; k < kStrain_; ++k) {

        tmpVec.push_back(rBernoulli(currentPlaf));

      }

      currentHap_.push_back(tmpVec);

    }

  }

  assert(currentHap_.size() == dEploidIO_->getPlaf().size());

}


double kgd::MCMCBASE::rBernoulli(double p) {

  double u = random_unit_.generate(entropy_source_.generator());

  return (u < p) ? 1.0 : 0.0;

}


void kgd::MCMCBASE::initializeExpectedWsaf(const std::vector<double> &proportion) {

  assert(currentExpectedWsaf_.size() == 0);

  currentExpectedWsaf_ = calcExpectedWsaf(proportion);

  assert(currentExpectedWsaf_.size() == nLoci_);

  cumExpectedWsaf_ = currentExpectedWsaf_;

}


std::vector<double> kgd::MCMCBASE::calcExpectedWsaf(const std::vector<double> &proportion) {
  //assert ( sumOfVec(proportion) == 1.0); // this fails ...

  std::vector<double> expectedWsaf;

  for (auto allele_site : currentHap_) {

    assert(kStrain_ == allele_site.size());

    double expectedWsaf_i = 0.0;

    for (size_t k = 0; k < kStrain(); k++) {

      expectedWsaf_i += allele_site[k] * proportion[k];

    }

    assert (expectedWsaf_i >= 0);

    expectedWsaf.push_back(expectedWsaf_i);

  }

  return expectedWsaf;

}


void kgd::MCMCBASE::computeDiagnostics() {

  dEploidIO_->setacceptRatio(acceptCount() / static_cast<double>(total_MCMC_iterations()));

  // average cumulate expectedWSAF
  for (size_t i = 0; i < cumExpectedWsaf_.size(); ++i) {

    cumExpectedWsaf_[i] /= static_cast<double>(MCMCParameters_.McmcSample());

  }

  std::vector<double> tmpLLKs1 = Utility::calcLLKs(dEploidIO_->getRefCount(),
                                                   dEploidIO_->getAltCount(),
                                                   cumExpectedWsaf_,
                                                   0,
                                                   cumExpectedWsaf_.size(),
                                                   MCMCParameters_.proposalUpdateScaling());

  dEploidIO_->setmeanThetallks(Utility::sumOfVec(tmpLLKs1));

  std::vector<double> wsaf_vec;

  for (size_t i = 0; i < nLoci(); ++i) {

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
                                                  MCMCParameters_.proposalUpdateScaling());

  dEploidIO_->setmaxLLKs(Utility::sumOfVec(tmpLLKs));

  double sum = std::accumulate(mcmcSample_->getSumLLKs().begin(), mcmcSample_->getSumLLKs().end(), 0.0);

  double mean = sum / mcmcSample_->getSumLLKs().size();

  double sq_sum = std::inner_product(mcmcSample_->getSumLLKs().begin(), mcmcSample_->getSumLLKs().end(),
                                     mcmcSample_->getSumLLKs().begin(), 0.0);

  double varLLKs = sq_sum / mcmcSample_->getSumLLKs().size() - mean * mean;
  double stdev = std::sqrt(varLLKs);

  dEploidIO_->setmeanllks(mean);
  dEploidIO_->setstdvllks(stdev);

  double dicByVar = (-2 * mean) + 4 * varLLKs / 2;
  dEploidIO_->setdicByVar(dicByVar);

  double dicWSAFBar = -2 * Utility::sumOfVec(tmpLLKs1);
  double dicByTheta = (-2 * mean) + (-2 * mean) - dicWSAFBar;

  dEploidIO_->setdicByTheta(dicByTheta);


}

