//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_subsetupdater.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::SubsetRelRateUpdater::SubsetRelRateUpdater(std::shared_ptr<Model> model_ptr) {
  //std::cout << "Creating a SubsetRelRateUpdater" << std::endl;
  DirichletUpdater::clear();
  name("Subset Relative Rates");
  _model = model_ptr;

}


kpl::SubsetRelRateUpdater::~SubsetRelRateUpdater() {
  //std::cout << "Destroying a SubsetRelRateUpdater" << std::endl;
}


double kpl::SubsetRelRateUpdater::calcLogPrior() {

  const Model::subset_sizes_t & subset_sizes = _model->getSubsetSizes();
  double log_num_sites = std::log(_model->getNumSites());
  unsigned num_subsets = _model->getNumSubsets();
  double log_prior = DirichletUpdater::calcLogPrior();

  for (unsigned i = 0; i < num_subsets-1; i++) {

    log_prior += std::log(subset_sizes[i]) - log_num_sites;

  }

  return log_prior;

}


void kpl::SubsetRelRateUpdater::pullFromModel() {

  const Model::subset_relrate_vect_t& relative_rates = _model->getSubsetRelRates();
  const Model::subset_sizes_t& subset_sizes   = _model->getSubsetSizes();

  unsigned num_sites   = _model->getNumSites();
  unsigned num_subsets = _model->getNumSubsets();

  assert(subset_sizes.size() == num_subsets);
  assert(relative_rates.size() == num_subsets);

  _curr_point.resize(num_subsets);

  for (unsigned i = 0; i < num_subsets; i++) {

    _curr_point[i] = relative_rates[i]*subset_sizes[i]/num_sites;

  }

}


void kpl::SubsetRelRateUpdater::pushToModel() {

  const Model::subset_sizes_t& subset_sizes = _model->getSubsetSizes();
  unsigned num_sites   = _model->getNumSites();
  unsigned num_subsets = _model->getNumSubsets();
  point_t tmp(num_subsets);

  for (unsigned i = 0; i < num_subsets; i++) {

    tmp[i] = _curr_point[i]*num_sites/subset_sizes[i];

  }

  _model->setSubsetRelRates(tmp, false);

}