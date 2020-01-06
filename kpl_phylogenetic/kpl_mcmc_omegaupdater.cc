//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_omegaupdater.h"


namespace kpl = kellerberrin::phylogenetic;



// member function bodies go here

kpl::OmegaUpdater::OmegaUpdater(QMatrix::SharedPtr q) {
  //std::cout << "OmegaUpdater being created" << std::endl;
  clear();
  name("Omega");
  assert(q);
  _q = q;

}


kpl::OmegaUpdater::~OmegaUpdater() {
  //std::cout << "OmegaUpdater being destroyed" << std::endl;
  _q.reset();

}


void kpl::OmegaUpdater::clear() {

  Updater::clear();
  _prev_point = 0.0;
  _q = nullptr;
  reset();

}


double kpl::OmegaUpdater::getCurrentPoint() const {

  return *(_q->getOmegaSharedPtr());

}


double kpl::OmegaUpdater::calcLogPrior() {

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


void kpl::OmegaUpdater::revert() {

  _q->setOmega(_prev_point);

}


void kpl::OmegaUpdater::proposeNewState() {
  // Save copy of _curr_point in case revert is necessary.
  _prev_point = getCurrentPoint();

  // Propose new value using multiplier with boldness _lambda
  double m = exp(lambda() * (lot()->uniform() - 0.5));

  _q->setOmega(m*_prev_point);

  // Calculate log of Hastings ratio
  logHastingsRatio(log(m));

  // This proposal invalidates all transition matrices and partials
  treeManipulator()->selectAllPartials();
  treeManipulator()->selectAllTMatrices();

}

