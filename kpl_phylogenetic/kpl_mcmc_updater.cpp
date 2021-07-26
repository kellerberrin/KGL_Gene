//
// Created by kellerberrin on 16/12/19.
//

#include "kel_exec_env.h"
#include "kpl_mcmc_updater.h"


namespace kpl = kellerberrin::phylogenetic;




kpl::Updater::Updater() {

  //std::cout << "Updater constructor called" << std::endl;
  clear();

}



void kpl::Updater::clear() {

  name_                   = "updater";
  tuning_                 = true;
  lambda_                 = 0.0001;
  weight_                 = 1.0;
  prob_                   = 0.0;
  target_acceptance_      = 0.3;
  n_accepts_               = 0;
  n_attempts_              = 0;
  heating_power_          = 1.0;
  prior_parameters_.clear();
  reset();

}


void kpl::Updater::reset() {

  log_hastings_ratio_ = 0.0;
  logJacobian(0.0);

}


void kpl::Updater::setLikelihood(Likelihood::SharedPtr likelihood) {

  likelihood_ptr_ = likelihood;

}


void kpl::Updater::setTreeManip(TreeManip::SharedPtr treemanip) {

  tree_manipulator_ptr_ = treemanip;

}




void kpl::Updater::setLot(Lot::SharedPtr lot) {

  lot_ptr_ = lot;

}


void kpl::Updater::setHeatingPower(double p) {

  heating_power_ = p;

}


void kpl::Updater::setLambda(double lambda) {

  lambda_ = lambda;

}


void kpl::Updater::setTuning(bool do_tune) {

  tuning_ = do_tune;
  n_accepts_ = 0;
  n_attempts_ = 0;

}


void kpl::Updater::tune(bool accepted) {

  n_attempts_++;

  if (tuning_) {

    double gamma_n = 10.0/(100.0 + (double)n_attempts_);

    if (accepted) {

      lambda_ *= 1.0 + gamma_n * (1.0 - target_acceptance_) / (2.0 * target_acceptance_);

    }
    else {

      lambda_ *= 1.0 - gamma_n * 0.5;

    }

    // Prevent run-away increases in boldness for low-information marginal densities
    if (lambda_ > 1000.0) {

      lambda_ = 1000.0;

    }

  }

}


void kpl::Updater::setTargetAcceptanceRate(double target) {

  target_acceptance_ = target;

}


void kpl::Updater::setPriorParameters(const std::vector<double> & c) {

  prior_parameters_.clear();
  prior_parameters_.assign(c.begin(), c.end());

}


void kpl::Updater::setWeight(double w) {

  weight_ = w;

}


void kpl::Updater::calcProb(double wsum) {

  assert(wsum > 0.0);
  prob_ = weight_ / wsum;

}


double kpl::Updater::getLambda() const {

  return lambda_;

}



double kpl::Updater::getWeight() const {

  return weight_;

}


double kpl::Updater::getAcceptPct() const {

  return (n_attempts_ == 0 ? 0.0 : (100.0 * n_accepts_ / n_attempts_));

}


double kpl::Updater::getNumUpdates() const {

  return n_attempts_;

}


std::string kpl::Updater::getUpdaterName() const {

  return name_;

}


double kpl::Updater::calcLogLikelihood() const {

  return likelihood_ptr_->calcLogLikelihood(*(tree_manipulator_ptr_->getTree()));

}


double kpl::Updater::update(double prev_lnL) {


  double prev_log_prior = calcLogPrior();

  // Clear any nodes previously selected so that we can detect those nodes
  // whose partials and/or transition_ probabilities need to be recalculated
  tree_manipulator_ptr_->deselectAllPartials();
  tree_manipulator_ptr_->deselectAllTMatrices();

  // Set model to proposed FSM_State and calculate log_hastings_ratio_
  proposeNewState();

  // Use alternative partials and transition_ probability buffer for any selected nodes
  // This allows us to easily revert to the previous values if the move is rejected
  tree_manipulator_ptr_->flipPartialsAndTMatrices();

  // Calculate the log-likelihood and log-prior for the proposed FSM_State
  double log_likelihood = calcLogLikelihood();
  double log_prior = calcLogPrior();

  // Decide whether to accept or reject the proposed FSM_State
  bool accept = true;
  double log_R = 0.0;
  double logu = 0.0;

  if (log_prior > LOG_ZERO_) {

    log_R = 0.0;
    log_R += heating_power_ * (log_likelihood - prev_lnL);
    log_R += heating_power_ * (log_prior - prev_log_prior);
    log_R += log_hastings_ratio_;
    log_R += logJacobian();

    logu = lot_ptr_->logUniform();

    if (logu > log_R) {

      accept = false;

    }

  }
  else {

    accept = false;

  }

  if (accept) {

    n_accepts_++;

  }
  else {

    revert();
    tree_manipulator_ptr_->flipPartialsAndTMatrices();
    log_likelihood = prev_lnL;

  }

  tune(accept);
  reset();


  return log_likelihood;

}


double kpl::Updater::calcEdgeLengthPrior() const {

  std::shared_ptr<Tree> tree = tree_manipulator_ptr_->getTree();
  assert(tree);

  double TL = tree_manipulator_ptr_->calcTreeLength();
  double num_edges = tree_manipulator_ptr_->countEdges();

  assert(prior_parameters_.size() == 3);
  double a = prior_parameters_[0];    // shape of Gamma prior on TL
  double b = prior_parameters_[1];    // scale of Gamma prior on TL
  double c = prior_parameters_[2];    // parameter of Dirichlet prior on edge length proportions

  // Calculate Gamma prior on tree length (TL)
  double log_gamma_prior_on_TL = (a - 1.0)*log(TL) - TL/b - a*log(b) - std::lgamma(a);

  // Calculate Dirichlet prior on edge length proportions
  //
  // Note that, for n edges, the Dirichlet prior density is
  //
  // p1^{c-1} p2^{c-1} ... pn^{c-1}
  // ------------------------------
  //    n*Gamma(c) / Gamma(n*c)
  //
  // where n = num_edges, pk = edge length k / TL and Gamma is the Gamma function.
  // If c == 1, then both numerator and denominator equal 1, so it is pointless
  // do loop over edge lengths.
  double log_edge_length_proportions_prior = std::lgamma(num_edges*c);

  if (c != 1.0) {

    for (auto nd : tree->getConstPreOrder()) {

      double edge_length_proportion = nd->getEdgeLength()/TL;
      log_edge_length_proportions_prior += (c - 1.0)*log(edge_length_proportion);

    }

    log_edge_length_proportions_prior -= std::lgamma(c)*num_edges;

  }

  double log_prior = log_gamma_prior_on_TL + log_edge_length_proportions_prior;

  return log_prior;

}


double kpl::Updater::getLogZero() {

  return LOG_ZERO_;

}


double kpl::Updater::calcLogTopologyPrior() const {

  std::shared_ptr<Tree> tree = tree_manipulator_ptr_->getTree();
  assert(tree);

  if (tree->isRooted()) {

    topo_prior_calculator_.chooseRooted();

  }
  else {

    topo_prior_calculator_.chooseUnrooted();

  }

  topo_prior_calculator_.setNTax(tree->numLeaves());

  unsigned m = tree->numInternals();

  double log_topology_prior = topo_prior_calculator_.getLogNormalizedTopologyPrior(m);

  return log_topology_prior;

}


void kpl::Updater::setTopologyPriorOptions(bool resclass, double C) {

  topo_prior_calculator_.setC(C);

  if (resclass) {

    topo_prior_calculator_.chooseResolutionClassPrior();

  }
  else {

    topo_prior_calculator_.choosePolytomyPrior();

  }

}