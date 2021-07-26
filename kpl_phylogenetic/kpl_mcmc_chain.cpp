//
// Created by kellerberrin on 16/12/19.
//


#include "kpl_mcmc_chain.h"
#include "kpl_mcmc_exchangeupdater.h"
#include "kpl_mcmc_omegaupdater.h"
#include "kpl_mcmc_pinvarupdater.h"
#include "kpl_mcmc_subsetupdater.h"
#include "kpl_mcmc_statefrequpdater.h"
#include "kel_exec_env.h"


namespace kpl = kellerberrin::phylogenetic;



kpl::Chain::Chain() {

  startTuning();

}



void kpl::Chain::startTuning() {

  for (auto& u : updaters_) {

    u->setTuning(true);

  }

}


void kpl::Chain::stopTuning() {

  for (auto u : updaters_) {

    u->setTuning(false);

  }
}


void kpl::Chain::setTreeFromNewick(const std::string& newick) {

  assert(updaters_.size() > 0);

  if (!tree_manip_ptr_) {

    tree_manip_ptr_.reset(new TreeManip);

  }

  tree_manip_ptr_->buildFromNewick(newick, false, true);

  for (auto u : updaters_) {

    u->setTreeManip(tree_manip_ptr_);

  }

}


size_t kpl::Chain::createUpdaters(const std::shared_ptr<Model>& model_ptr,
                                  const std::shared_ptr<Lot>& lot,
                                  const std::shared_ptr<Likelihood>& likelihood) {

  model_ptr_ = model_ptr;
  lot_ptr_ = lot;
  updaters_.clear();

  double wstd             = 1.0;
  double sum_weights      = 0.0;
  double wtreelength      = 1.0;
  double wtreetopology    = 19.0;
  double wpolytomy        = 0.0;

  if (model_ptr_->isAllowPolytomies()) {

    wstd             = 1.0;
    wtreelength      = 2.0;
    wtreetopology    = 10.0;
    wpolytomy        = 5.0;

  }

  // Add FSM_State frequency parameter updaters to updaters_
  const Model::state_freq_params_t& statefreq_shptr_vect = model_ptr_->getStateFreqParams();
  for (auto& statefreq_shptr : statefreq_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<StateFreqUpdater>(statefreq_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(0.001);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters(std::vector<double>(statefreq_shptr->getStateFreqsSharedPtr()->size(), 1.0));
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    updaters_.push_back(local_updater_ptr);

  }

  // Add exchangeability parameter updaters to updaters_
  const Model::exchangeability_params_t& exchangeability_shptr_vect = model_ptr_->getExchangeabilityParams();
  for (auto& exchangeability_shptr : exchangeability_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<ExchangeabilityUpdater>(exchangeability_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    updaters_.push_back(local_updater_ptr);

  }

  // Add rate variance parameter updaters to updaters_
  const Model::ratevar_params_t & ratevar_shptr_vect = model_ptr_->getRateVarParams();
  for (auto& ratevar_shptr : ratevar_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<GammaRateVarUpdater>(ratevar_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    updaters_.push_back(local_updater_ptr);

  }

  // Add pinvar parameter updaters to updaters_
  const Model::pinvar_params_t& pinvar_shptr_vect = model_ptr_->getPinvarParams();
  for (auto& pinvar_shptr : pinvar_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<PinvarUpdater>(pinvar_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(0.5);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    updaters_.push_back(local_updater_ptr);

  }

  // Add omega parameter updaters to updaters_
  const Model::omega_params_t& omega_shptr_vect = model_ptr_->getOmegaParams();
  std::cout << "Num Omega Params:" << omega_shptr_vect.size() << std::endl;
  for (auto& omega_shptr : omega_shptr_vect) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<OmegaUpdater>(omega_shptr);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters({1.0, 1.0});
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    updaters_.push_back(local_updater_ptr);

  }

  // Add subset relative rate parameter updater to updaters_
  if (not model_ptr_->isFixedSubsetRelRates()) {

    Updater::SharedPtr local_updater_ptr = std::make_shared<SubsetRelRateUpdater>(model_ptr_);
    local_updater_ptr->setLikelihood(likelihood);
    local_updater_ptr->setLot(lot);
    local_updater_ptr->setLambda(1.0);
    local_updater_ptr->setTargetAcceptanceRate(0.3);
    local_updater_ptr->setPriorParameters(std::vector<double>(model_ptr_->getNumSubsets(), 1.0));
    local_updater_ptr->setWeight(wstd); sum_weights += wstd;
    updaters_.push_back(local_updater_ptr);

  }

  // Add tree updater and tree length updater to updaters_
  if (!model_ptr_->isFixedTree()) {

    double tree_length_shape = 1.0;
    double tree_length_scale = 10.0;
    double dirichlet_param   = 1.0;

    Updater::SharedPtr updater_ptr = std::make_shared<TreeUpdater>();
    updater_ptr->setLikelihood(likelihood);
    updater_ptr->setLot(lot);
    updater_ptr->setLambda(0.5);
    updater_ptr->setTargetAcceptanceRate(0.3);
    updater_ptr->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
    updater_ptr->setTopologyPriorOptions(model_ptr_->isResolutionClassTopologyPrior(), model_ptr_->getTopologyPriorC());
    updater_ptr->setWeight(wtreetopology); sum_weights += wtreetopology;

    updaters_.push_back(updater_ptr);

    if (model_ptr_->isAllowPolytomies()) {

      Updater::SharedPtr local_updater_ptr = std::make_shared<PolytomyUpdater>();
      local_updater_ptr->setLikelihood(likelihood);
      local_updater_ptr->setLot(lot);
      local_updater_ptr->setLambda(0.05);
      local_updater_ptr->setTargetAcceptanceRate(0.3);
      local_updater_ptr->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
      local_updater_ptr->setTopologyPriorOptions(model_ptr_->isResolutionClassTopologyPrior(), model_ptr_->getTopologyPriorC());
      local_updater_ptr->setWeight(wpolytomy); sum_weights += wpolytomy;

      updaters_.push_back(local_updater_ptr);

    }

    updater_ptr = std::make_shared<TreeLengthUpdater>();
    updater_ptr->setLikelihood(likelihood);
    updater_ptr->setLot(lot);
    updater_ptr->setLambda(0.2);
    updater_ptr->setTargetAcceptanceRate(0.3);
    updater_ptr->setPriorParameters({tree_length_shape, tree_length_scale, dirichlet_param});
    updater_ptr->setTopologyPriorOptions(model_ptr_->isResolutionClassTopologyPrior(), model_ptr_->getTopologyPriorC());
    updater_ptr->setWeight(wtreelength); sum_weights += wtreelength;

    updaters_.push_back(updater_ptr);

  }

  for (auto& updater : updaters_) {

    updater->calcProb(sum_weights);

  }

  return updaters_.size();

}


void kpl::Chain::setHeatingPower(double p) {

  heating_power_ = p;
  for (auto u : updaters_) {

    u->setHeatingPower(p);

  }

}


double kpl::Chain::getChainIndex() const {

  return chain_index_;

}


void kpl::Chain::setChainIndex(unsigned idx) {

  chain_index_ = idx;

}


std::vector<std::string> kpl::Chain::getUpdaterNames() const {

  std::vector<std::string> v;

  for (auto u : updaters_) {

    v.push_back(u->getUpdaterName());

  }

  return v;

}


std::vector<double> kpl::Chain::getAcceptPercentages() const {

  std::vector<double> v;

  for (auto u : updaters_) {

    v.push_back(u->getAcceptPct());

  }

  return v;

}


std::vector<unsigned> kpl::Chain::getNumUpdates() const {

  std::vector<unsigned> v;

  for (auto u : updaters_) {

    v.push_back(u->getNumUpdates());

  }

  return v;

}


std::vector<double> kpl::Chain::getLambdas() const {

  std::vector<double> v;

  for (auto u : updaters_) {

    v.push_back(u->getLambda());

  }

  return v;

}


void kpl::Chain::setLambdas(std::vector<double> & v) {

  assert(v.size() == updaters_.size());

  unsigned index = 0;
  for (auto u : updaters_) {

    u->setLambda(v[index++]);

  }

}


double kpl::Chain::calcLogLikelihood() const {

  return updaters_[0]->calcLogLikelihood();

}


double kpl::Chain::calcLogJointPrior() const {

  double lnP = 0.0;

  for (auto u : updaters_) {

    if (u->name() != "Tree Length" && u->name() != "Polytomies" ) {

      lnP += u->calcLogPrior();

    }

  }

  return lnP;

}


void kpl::Chain::start() {

  tree_manip_ptr_->selectAllPartials();
  tree_manip_ptr_->selectAllTMatrices();
  log_likelihood_ = calcLogLikelihood();
  tree_manip_ptr_->deselectAllPartials();
  tree_manip_ptr_->deselectAllTMatrices();

}

void kpl::Chain::stop() {

}


void kpl::Chain::nextStep(int iteration) {

  double u = lot_ptr_->uniform();
  double cumprob{0.0};
  size_t i{0};

  for (auto const& updater : updaters_) {

    cumprob += updater->probability();

    if (u <= cumprob) {

      break;

    }

    i++;

  }

  if (i >= updaters_.size()) {

    ExecEnv::log().error("Chain::nextStep; updater index: {} is invalid, updater count: {}", i, updaters_.size());
    i = 0;

  }

  log_likelihood_ = updaters_[i]->update(log_likelihood_);

}


