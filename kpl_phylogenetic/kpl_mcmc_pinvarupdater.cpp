//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_pinvarupdater.h"


namespace kpl = kellerberrin::phylogenetic;



// member function bodies go here

kpl::PinvarUpdater::PinvarUpdater(std::shared_ptr<ASRV> asrv) {

  clear();
  name("Proportion of Invariable Sites");
  assert(asrv);
  _asrv = asrv;

}


kpl::PinvarUpdater::~PinvarUpdater() {

  _asrv.reset();

}


void kpl::PinvarUpdater::clear() {

  Updater::clear();
  _prev_point = 0.0;
  _asrv = nullptr;
  reset();

}


double kpl::PinvarUpdater::getCurrentPoint() const {

  return *(_asrv->getPinvarSharedPtr());

}


double kpl::PinvarUpdater::calcLogPrior() {
  // Assumes Beta(a,b) prior with mean a/(a+b) and variance a*b/((a + b + 1)*(a + b)^2)
  assert(priorParameters().size() == 2);
  double prior_a = priorParameters()[0];
  double prior_b = priorParameters()[1];

  double log_prior = 0.0;
  double curr_point = getCurrentPoint();

  if (curr_point > 0.0 && curr_point < 1.0) {

    log_prior += (prior_a - 1.0)*std::log(curr_point);
    log_prior += (prior_b - 1.0)*std::log(1.0 - curr_point);
    log_prior += std::lgamma(prior_a + prior_b);
    log_prior -= std::lgamma(prior_a);
    log_prior -= std::lgamma(prior_b);

  }
  else {

    log_prior = Updater::logZero();

  }

  return log_prior;

}


void kpl::PinvarUpdater::revert() {

  _asrv->setPinvar(_prev_point);

}


void kpl::PinvarUpdater::proposeNewState() {

  if (lambda() > 1.0) {

    lambda(1.0);

  }

  // Save copy of _curr_point in case revert is necessary.
  _prev_point = getCurrentPoint();

  // Propose new value using uniform window of width _lambda centered over _prev_point
  double p = (_prev_point - lambda() /2.0) + lambda() * lot()->uniform();

  if (p < 0.0) {

    p = -p;

  }
  else if (p > 1.0) {

    p = 1.0 - (p - 1.0);

  }

  _asrv->setPinvar(p);

  logHastingsRatio(0.0);  //symmetric proposal

  // This proposal invalidates all transition_ matrices and partials
  treeManipulator()->selectAllPartials();
  treeManipulator()->selectAllTMatrices();

}

