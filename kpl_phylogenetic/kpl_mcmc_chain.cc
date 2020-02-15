//
// Created by kellerberrin on 16/12/19.
//

#include "kpl_mcmc_chain.h"
#include "kpl_mcmc_exchangeupdater.h"
#include "kpl_mcmc_omegaupdater.h"
#include "kpl_mcmc_pinvarupdater.h"
#include "kpl_mcmc_subsetupdater.h"
#include "kpl_mcmc_statefrequpdater.h"



namespace kpl = kellerberrin::phylogenetic;



kpl::Chain::Chain() {
  //std::cout << "Chain being created" << std::endl;
  clear();
}


kpl::Chain::~Chain() {
  //std::cout << "Chain being destroyed" << std::endl;
}


void kpl::Chain::clear() {

  _log_likelihood = 0.0;
  _updaters.clear();
  _chain_index = 0;
  setHeatingPower(1.0);
  startTuning();

}


void kpl::Chain::startTuning() {

  for (auto u : _updaters) {

    u->setTuning(true);

  }

}


void kpl::Chain::stopTuning() {

  for (auto u : _updaters) {

    u->setTuning(false);

  }
}


void kpl::Chain::setTreeFromNewick(std::string & newick) {

  assert(_updaters.size() > 0);

  if (!_tree_manipulator) {

    _tree_manipulator.reset(new TreeManip);

  }

  _tree_manipulator->buildFromNewick(newick, false, true);

  for (auto u : _updaters) {

    u->setTreeManip(_tree_manipulator);

  }

}


unsigned kpl::Chain::createUpdaters(std::shared_ptr<Model> model_ptr, Lot::SharedPtr lot, Likelihood::SharedPtr likelihood) {

  _model = model_ptr;
  _lot = lot;
  _updaters.clear();

  double wstd             = 1.0;
  double sum_weights      = 0.0;
  double wtreelength      = 1.0;
  double wtreetopology    = 19.0;
  double wpolytomy        = 0.0;

  if (_model->isAllowPolytomies()) {

    wstd             = 1.0;
    wtreelength      = 2.0;
    wtreetopology    = 10.0;
    wpolytomy        = 5.0;

  }

  // Add state frequency parameter updaters to _updaters
  const Model::state_freq_params_t& statefreq_shptr_vect = _model->getStateFreqParams();
  for (auto& statefreq_shptr : statefreq_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<StateFreqUpdater>(statefreq_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(0.001);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters(std::vector<double>(statefreq_shptr->getStateFreqsSharedPtr()->size(), 1.0));
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    _updaters.push_back(local_updater_ptr);

  }

  // Add exchangeability parameter updaters to _updaters
  const Model::exchangeability_params_t& exchangeability_shptr_vect = _model->getExchangeabilityParams();
  for (auto& exchangeability_shptr : exchangeability_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<ExchangeabilityUpdater>(exchangeability_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    _updaters.push_back(local_updater_ptr);

  }

  // Add rate variance parameter updaters to _updaters
  const Model::ratevar_params_t & ratevar_shptr_vect = _model->getRateVarParams();
  for (auto& ratevar_shptr : ratevar_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<GammaRateVarUpdater>(ratevar_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    _updaters.push_back(local_updater_ptr);

  }

  // Add pinvar parameter updaters to _updaters
  const Model::pinvar_params_t& pinvar_shptr_vect = _model->getPinvarParams();
  for (auto& pinvar_shptr : pinvar_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<PinvarUpdater>(pinvar_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(0.5);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    _updaters.push_back(local_updater_ptr);

  }

  // Add omega parameter updaters to _updaters
  const Model::omega_params_t& omega_shptr_vect = _model->getOmegaParams();
  std::cout << "Num Omega Params:" << omega_shptr_vect.size() << std::endl;
  for (auto& omega_shptr : omega_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<OmegaUpdater>(omega_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    _updaters.push_back(local_updater_ptr);

  }

  // Add subset relative rate parameter updater to _updaters
  if (!_model->isFixedSubsetRelRates()) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<SubsetRelRateUpdater>(_model);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters(std::vector<double>(_model->getNumSubsets(), 1.0));
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    _updaters.push_back(local_updater_ptr);

  }

  // Add tree updater and tree length updater to _updaters
  if (!_model->isFixedTree()) {

    double tree_length_shape = 1.0;
    double tree_length_scale = 10.0;
    double dirichlet_param   = 1.0;

    Updater::SharedPtr updater_ptr = std::make_shared<TreeUpdater>();
    updater_ptr->setLikelihood(likelihood);
    updater_ptr->setLot(lot);
    updater_ptr->setLambda(0.5);
    updater_ptr->setTargetAcceptanceRate(0.3);
    updater_ptr->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
    updater_ptr->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
    updater_ptr->setWeight(wtreetopology); sum_weights += wtreetopology;

    _updaters.push_back(updater_ptr);

    if (_model->isAllowPolytomies()) {

      Updater::SharedPtr local_updater_ptr = std::make_shared<PolytomyUpdater>();
      local_updater_ptr->setLikelihood(likelihood);
      local_updater_ptr->setLot(lot);
      local_updater_ptr->setLambda(0.05);
      local_updater_ptr->setTargetAcceptanceRate(0.3);
      local_updater_ptr->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
      local_updater_ptr->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
      local_updater_ptr->setWeight(wpolytomy); sum_weights += wpolytomy;

      _updaters.push_back(local_updater_ptr);

    }

    updater_ptr = std::make_shared<TreeLengthUpdater>();
    updater_ptr->setLikelihood(likelihood);
    updater_ptr->setLot(lot);
    updater_ptr->setLambda(0.2);
    updater_ptr->setTargetAcceptanceRate(0.3);
    updater_ptr->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
    updater_ptr->setTopologyPriorOptions(_model->isResolutionClassTopologyPrior(), _model->getTopologyPriorC());
    updater_ptr->setWeight(wtreelength); sum_weights += wtreelength;

    _updaters.push_back(updater_ptr);

  }

  for (auto& updater : _updaters) {

    updater->calcProb(sum_weights);

  }

  return (unsigned)_updaters.size();

}



kpl::TreeManip::SharedPtr kpl::Chain::getTreeManip() {

  return _tree_manipulator;

}


std::shared_ptr<kpl::Model> kpl::Chain::getModel() {

  return _model;

}


double kpl::Chain::getLogLikelihood() const {

  return _log_likelihood;

}


double kpl::Chain::getHeatingPower() const {

  return _heating_power;

}


void kpl::Chain::setHeatingPower(double p) {

  _heating_power = p;
  for (auto u : _updaters) {

    u->setHeatingPower(p);

  }

}


double kpl::Chain::getChainIndex() const {

  return _chain_index;

}


void kpl::Chain::setChainIndex(unsigned idx) {

  _chain_index = idx;

}


kpl::Updater::SharedPtr kpl::Chain::findUpdaterByName(std::string name) {

  Updater::SharedPtr retval = nullptr;

  for (auto u : _updaters) {

    if (u->getUpdaterName() == name) {

      retval = u;
      break;

    }

  }

  assert(retval != nullptr);
  return retval;

}


std::vector<std::string> kpl::Chain::getUpdaterNames() const {

  std::vector<std::string> v;

  for (auto u : _updaters) {

    v.push_back(u->getUpdaterName());

  }

  return v;

}


std::vector<double> kpl::Chain::getAcceptPercentages() const {

  std::vector<double> v;

  for (auto u : _updaters) {

    v.push_back(u->getAcceptPct());

  }

  return v;

}


std::vector<unsigned> kpl::Chain::getNumUpdates() const {

  std::vector<unsigned> v;

  for (auto u : _updaters) {

    v.push_back(u->getNumUpdates());

  }

  return v;

}


std::vector<double> kpl::Chain::getLambdas() const {

  std::vector<double> v;

  for (auto u : _updaters) {

    v.push_back(u->getLambda());

  }

  return v;

}


void kpl::Chain::setLambdas(std::vector<double> & v) {

  assert(v.size() == _updaters.size());

  unsigned index = 0;
  for (auto u : _updaters) {

    u->setLambda(v[index++]);

  }

}


double kpl::Chain::calcLogLikelihood() const {

  return _updaters[0]->calcLogLikelihood();

}


double kpl::Chain::calcLogJointPrior() const {

  double lnP = 0.0;

  for (auto u : _updaters) {

    if (u->name() != "Tree Length" && u->name() != "Polytomies" ) {

      lnP += u->calcLogPrior();

    }

  }

  return lnP;

}


void kpl::Chain::start() {

  _tree_manipulator->selectAllPartials();
  _tree_manipulator->selectAllTMatrices();
  _log_likelihood = calcLogLikelihood();
  _tree_manipulator->deselectAllPartials();
  _tree_manipulator->deselectAllTMatrices();

}

void kpl::Chain::stop() {

}


void kpl::Chain::nextStep(int iteration) {

  assert(_lot);

  double u = _lot->uniform();
  double cumprob = 0.0;
  unsigned i = 0;

  for (auto updater : _updaters) {

    cumprob += updater->probability();

    if (u <= cumprob) {

      break;

    }

    i++;

  }

  assert(i < _updaters.size());

  _log_likelihood = _updaters[i]->update(_log_likelihood);

}


