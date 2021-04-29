//
// Created by kellerberrin on 16/12/19.
//

#include "kpl_mcmc_gamma.h"


namespace kpl = kellerberrin::phylogenetic;



kpl::GammaRateVarUpdater::GammaRateVarUpdater(std::shared_ptr<ASRV> asrv) {

  //std::cout << "GammaRateVarUpdater being created" << std::endl;
  clear();
  name("Gamma Rate Variance");
  assert(asrv);
  _asrv = asrv;

}


kpl::GammaRateVarUpdater::~GammaRateVarUpdater() {
  //std::cout << "GammaRateVarUpdater being destroyed" << std::endl;

  _asrv.reset();

}


void kpl::GammaRateVarUpdater::clear() {

  Updater::clear();
  _prev_point = 0.0;
  _asrv = nullptr;
  reset();

}


double kpl::GammaRateVarUpdater::getCurrentPoint() const {

  return _asrv->getRateVar();

}


double kpl::GammaRateVarUpdater::calcLogPrior() {
  // Assumes Gamma(a,b) prior with mean a*b and variance a*b^2
  assert(priorParameters().size() == 2);

  double prior_a = priorParameters()[0];
  double prior_b = priorParameters()[1];

  double log_prior = 0.0;
  double curr_point = getCurrentPoint();

  if (curr_point > 0.0) {

    log_prior += (prior_a - 1.0)*std::log(curr_point);
    log_prior -= curr_point/prior_b;
    log_prior -= prior_a*std::log(prior_b);
    log_prior -= std::lgamma(prior_a);

  }
  else {

    log_prior = Updater::logZero();

  }

  return log_prior;

}


void kpl::GammaRateVarUpdater::revert() {

  _asrv->setRateVar(_prev_point);

}


void kpl::GammaRateVarUpdater::proposeNewState() {
  // Save copy of _curr_point in case revert is necessary.
  _prev_point = getCurrentPoint();

  // Propose new value using multiplier with boldness lambda_
  double m = exp(lambda() *(lot()->uniform() - 0.5));
  _asrv->setRateVar(m*_prev_point);

  // Calculate log of Hastings ratio
  logHastingsRatio(std::log(m));

  // This proposal invalidates all transition_ matrices and partials
  treeManipulator()->selectAllPartials();
  treeManipulator()->selectAllTMatrices();

}

