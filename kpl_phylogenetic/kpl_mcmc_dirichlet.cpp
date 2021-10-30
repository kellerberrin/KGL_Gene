//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_dirichlet.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::DirichletUpdater::DirichletUpdater() {
  // std::cout << "Creating DirichletUpdater object" << std::endl;
  clear();

}


kpl::DirichletUpdater::~DirichletUpdater() {
  // std::cout << "Destroying DirichletUpdater object" << std::endl;

}


void kpl::DirichletUpdater::clear() {

  Updater::clear();
  _prev_point.clear();

}


double kpl::DirichletUpdater::calcLogPrior() {

  pullFromModel();

  assert(_curr_point.size() > 0);
  assert(_curr_point.size() == priorParameters().size());

  bool flat_prior = true;
  bool bad_point = false;
  double log_prior = 0.0;
  double prior_param_sum = 0.0;

  for (unsigned i = 0; i < _curr_point.size(); ++i) {

    if (priorParameters()[i] != 1.0) {

      flat_prior = false;

    }
    if (_curr_point[i] == 0.0) {

      bad_point = true;

    }

    log_prior += (priorParameters()[i] - 1.0)*std::log(_curr_point[i]);
    log_prior -= std::lgamma(priorParameters()[i]);
    prior_param_sum += priorParameters()[i];

  }

  if (flat_prior) {

    return std::lgamma(prior_param_sum);

  }
  else if (bad_point) {

    return Updater::logZero();

  }
  else {

    log_prior += std::lgamma(prior_param_sum);

  }

  return log_prior;

}


void kpl::DirichletUpdater::proposeNewState() {
  // Save length of _curr_point.
  pullFromModel();
  unsigned dim = (unsigned)_curr_point.size();

  // Save copy of _curr_point in case revert is necessary.
  _prev_point.assign(_curr_point.begin(), _curr_point.end());

  // Determine parameters of Dirichlet forward proposal distribution and, at the same time,
  // draw gamma deviates that will be used to form the proposed point.
  std::vector<double> forward_params(dim, 0.0);
  for (unsigned i = 0; i < dim; ++i) {
    // Calculate ith forward parameter
    double alpha_i = 1.0 + _prev_point[i] / lambda();

    if (alpha_i < 1.e-12) {

      alpha_i = 1.e-12;

    }

    forward_params[i] = alpha_i;

    // Draw ith gamma deviate
    _curr_point[i] = 0.0;

    if (alpha_i > 0.0) {

      _curr_point[i] = lot()->gamma(alpha_i, 1.0);

    }

  }

  double sum_gamma_deviates     = std::accumulate(_curr_point.begin(), _curr_point.end(), 0.0);
  double sum_forward_parameters = std::accumulate(forward_params.begin(), forward_params.end(), 0.0);

  // Choose new FSM_State by sampling from forward proposal distribution.
  // We've already stored gamma deviates in _curr_point, now just need to normalize them.
  for (unsigned i = 0; i < dim; ++i) {

    _curr_point[i] /= sum_gamma_deviates;

  }

  // Determine probfailure density of the forward proposal
  double log_forward_density = 0.0;
  for (unsigned i = 0; i < dim; ++i) {

    log_forward_density += (forward_params[i] - 1.0) * std::log(_prev_point[i]);
    log_forward_density -= std::lgamma(forward_params[i]);

  }

  log_forward_density += std::lgamma(sum_forward_parameters);

  // Determine parameters of Dirichlet reverse proposal distribution
  std::vector<double> reverse_params(dim, 0.0);
  for (unsigned i = 0; i < dim; ++i) {

    reverse_params[i] = 1.0 + _curr_point[i] / lambda();

  }

  double sum_reverse_parameters = std::accumulate(reverse_params.begin(), reverse_params.end(), 0.0);

  // determine probfailure density of the reverse proposal
  double log_reverse_density = 0.0;
  for (unsigned i = 0; i < dim; ++i) {

    log_reverse_density += (reverse_params[i] - 1.0)*std::log(_curr_point[i]);
    log_reverse_density -= std::lgamma(reverse_params[i]);

  }
  log_reverse_density += std::lgamma(sum_reverse_parameters);

  // calculate the logarithm of the Hastings ratio
  logHastingsRatio(log_reverse_density - log_forward_density);

  pushToModel();

  // This proposal invalidates all transition_ matrices and partials
  treeManipulator()->selectAllPartials();
  treeManipulator()->selectAllTMatrices();

}


void kpl::DirichletUpdater::revert() {

  std::copy(_prev_point.begin(), _prev_point.end(), _curr_point.begin());
  pushToModel();

}