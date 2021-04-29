//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_subsetupdater.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::SubsetRelRateUpdater::SubsetRelRateUpdater(const std::shared_ptr<Model>& model_ptr) {

  name("Subset Relative Rates");
  model_ptr_ = model_ptr;

}



double kpl::SubsetRelRateUpdater::calcLogPrior() {

  const Model::subset_sizes_t & subset_sizes = model_ptr_->getSubsetSizes();
  double log_num_sites = std::log(model_ptr_->getNumSites());
  unsigned num_subsets = model_ptr_->getNumSubsets();
  double log_prior = DirichletUpdater::calcLogPrior();

  for (unsigned i = 0; i < num_subsets-1; i++) {

    log_prior += std::log(subset_sizes[i]) - log_num_sites;

  }

  return log_prior;

}


void kpl::SubsetRelRateUpdater::pullFromModel() {

  const Model::subset_relrate_vect_t& relative_rates = model_ptr_->getSubsetRelRates();
  const Model::subset_sizes_t& subset_sizes   = model_ptr_->getSubsetSizes();

  unsigned num_sites   = model_ptr_->getNumSites();
  unsigned num_subsets = model_ptr_->getNumSubsets();

  assert(subset_sizes.size() == num_subsets);
  assert(relative_rates.size() == num_subsets);

  _curr_point.resize(num_subsets);

  for (unsigned i = 0; i < num_subsets; i++) {

    _curr_point[i] = relative_rates[i]*subset_sizes[i]/num_sites;

  }

}


void kpl::SubsetRelRateUpdater::pushToModel() {

  const Model::subset_sizes_t& subset_sizes = model_ptr_->getSubsetSizes();
  unsigned num_sites   = model_ptr_->getNumSites();
  unsigned num_subsets = model_ptr_->getNumSubsets();
  point_t tmp(num_subsets);

  for (unsigned i = 0; i < num_subsets; i++) {

    tmp[i] = _curr_point[i]*num_sites/subset_sizes[i];

  }

  model_ptr_->setSubsetRelRates(tmp, false);

}