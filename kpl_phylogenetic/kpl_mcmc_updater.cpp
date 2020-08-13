//
// Created by kellerberrin on 16/12/19.
//

#include "kpl_mcmc_updater.h"


namespace kpl = kellerberrin::phylogenetic;




kpl::Updater::Updater() {

  //std::cout << "Updater constructor called" << std::endl;
  clear();

}


kpl::Updater::~Updater() {

  //std::cout << "Updater destructor called" << std::endl;

}


void kpl::Updater::clear() {

  _name                   = "updater";
  _tuning                 = true;
  _lambda                 = 0.0001;
  _weight                 = 1.0;
  _prob                   = 0.0;
  _target_acceptance      = 0.3;
  _naccepts               = 0;
  _nattempts              = 0;
  _heating_power          = 1.0;
  _prior_parameters.clear();
  reset();

}


void kpl::Updater::reset() {

  _log_hastings_ratio = 0.0;
  logJacobian(0.0);

}


void kpl::Updater::setLikelihood(Likelihood::SharedPtr likelihood) {

  _likelihood = likelihood;

}


void kpl::Updater::setTreeManip(TreeManip::SharedPtr treemanip) {

  _tree_manipulator = treemanip;

}


kpl::TreeManip::SharedPtr kpl::Updater::getTreeManip() const {

  return _tree_manipulator;

}


void kpl::Updater::setLot(Lot::SharedPtr lot) {

  _lot = lot;

}


void kpl::Updater::setHeatingPower(double p) {

  _heating_power = p;

}


void kpl::Updater::setLambda(double lambda) {

  _lambda = lambda;

}


void kpl::Updater::setTuning(bool do_tune) {

  _tuning = do_tune;
  _naccepts = 0;
  _nattempts = 0;

}


void kpl::Updater::tune(bool accepted) {

  _nattempts++;

  if (_tuning) {

    double gamma_n = 10.0/(100.0 + (double)_nattempts);

    if (accepted) {

      _lambda *= 1.0 + gamma_n*(1.0 - _target_acceptance)/(2.0*_target_acceptance);

    }
    else {

      _lambda *= 1.0 - gamma_n*0.5;

    }

    // Prevent run-away increases in boldness for low-information marginal densities
    if (_lambda > 1000.0) {

      _lambda = 1000.0;

    }

  }

}


void kpl::Updater::setTargetAcceptanceRate(double target) {

  _target_acceptance = target;

}


void kpl::Updater::setPriorParameters(const std::vector<double> & c) {

  _prior_parameters.clear();
  _prior_parameters.assign(c.begin(), c.end());

}


void kpl::Updater::setWeight(double w) {

  _weight = w;

}


void kpl::Updater::calcProb(double wsum) {

  assert(wsum > 0.0);
  _prob = _weight/wsum;

}


double kpl::Updater::getLambda() const {

  return _lambda;

}


double kpl::Updater::getProb() const {

  return _prob;

}


double kpl::Updater::getWeight() const {

  return _weight;

}


double kpl::Updater::getAcceptPct() const {

  return (_nattempts == 0 ? 0.0 : (100.0*_naccepts/_nattempts));

}


double kpl::Updater::getNumUpdates() const {

  return _nattempts;

}


std::string kpl::Updater::getUpdaterName() const {

  return _name;

}


double kpl::Updater::calcLogLikelihood() const {

  return _likelihood->calcLogLikelihood(_tree_manipulator->getTree());

}


double kpl::Updater::update(double prev_lnL) {


  double prev_log_prior = calcLogPrior();

  // Clear any nodes previously selected so that we can detect those nodes
  // whose partials and/or transition probabilities need to be recalculated
  _tree_manipulator->deselectAllPartials();
  _tree_manipulator->deselectAllTMatrices();

  // Set model to proposed state and calculate _log_hastings_ratio
  proposeNewState();

  // Use alternative partials and transition probability buffer for any selected nodes
  // This allows us to easily revert to the previous values if the move is rejected
  _tree_manipulator->flipPartialsAndTMatrices();

  // Calculate the log-likelihood and log-prior for the proposed state
  double log_likelihood = calcLogLikelihood();
  double log_prior = calcLogPrior();

  // Decide whether to accept or reject the proposed state
  bool accept = true;
  double log_R = 0.0;
  double logu = 0.0;

  if (log_prior > _log_zero) {

    log_R = 0.0;
    log_R += _heating_power*(log_likelihood - prev_lnL);
    log_R += _heating_power*(log_prior - prev_log_prior);
    log_R += _log_hastings_ratio;
    log_R += logJacobian();

    logu = _lot->logUniform();

    if (logu > log_R) {

      accept = false;

    }

  }
  else {

    accept = false;

  }

  if (accept) {

    _naccepts++;

  }
  else {

    revert();
    _tree_manipulator->flipPartialsAndTMatrices();
    log_likelihood = prev_lnL;

  }

  tune(accept);
  reset();


  return log_likelihood;

}


double kpl::Updater::calcEdgeLengthPrior() const {

  Tree::SharedPtr tree = _tree_manipulator->getTree();
  assert(tree);

  double TL = _tree_manipulator->calcTreeLength();
  double num_edges = _tree_manipulator->countEdges();

  assert(_prior_parameters.size() == 3);
  double a = _prior_parameters[0];    // shape of Gamma prior on TL
  double b = _prior_parameters[1];    // scale of Gamma prior on TL
  double c = _prior_parameters[2];    // parameter of Dirichlet prior on edge length proportions

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

  return _log_zero;

}


double kpl::Updater::calcLogTopologyPrior() const {

  Tree::SharedPtr tree = _tree_manipulator->getTree();
  assert(tree);

  if (tree->isRooted()) {

    _topo_prior_calculator.chooseRooted();

  }
  else {

    _topo_prior_calculator.chooseUnrooted();

  }

  _topo_prior_calculator.setNTax(tree->numLeaves());

  unsigned m = tree->numInternals();

  double log_topology_prior = _topo_prior_calculator.getLogNormalizedTopologyPrior(m);

  return log_topology_prior;

}


void kpl::Updater::setTopologyPriorOptions(bool resclass, double C) {

  _topo_prior_calculator.setC(C);

  if (resclass) {

    _topo_prior_calculator.chooseResolutionClassPrior();

  }
  else {

    _topo_prior_calculator.choosePolytomyPrior();

  }

}