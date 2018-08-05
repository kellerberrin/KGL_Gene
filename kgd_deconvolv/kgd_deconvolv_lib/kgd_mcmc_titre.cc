//
// Created by kellerberrin on 27/06/18.
//

#include <ctime>
#include <cstdio>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_mcmc_titre.h"
#include "kgd_deconvolv_app.h"


namespace kgd = kellerberrin::deconvolv;


kgd::MCMCTITRE::MCMCTITRE(size_t k_strains,
                          double mean_log_titre,
                          double sd_log_titre,
                          double update_scale) : k_strains_(k_strains) ,
                                                mean_log_titre_(mean_log_titre),
                                                sd_log_titre_(sd_log_titre),
                                                update_scale_(update_scale) {


  std_norm_ = std::make_shared<StdNormalDistribution>();
  entropy_source_ = std::make_shared<EntropySource>();

  randomizeTitre();

  ExecEnv::log().info("Titre initialization values: {}", titreText());
  ExecEnv::log().info("Proportion initialization values: {}", proportionsText());

}


void kgd::MCMCTITRE::randomizeProportions() {

  std::vector<double> proportions;

  for (size_t k = 0; k < kStrain(); ++k) {

    proportions.push_back(std::fabs(initialTitreNormalVariable()));

  }

  Utility::normalizeBySum(proportions);

  proportion2Titre(proportions);

}

void kgd::MCMCTITRE::randomizeTitre() {

  currentTitre_.clear();

  for (size_t k = 0; k < kStrain(); ++k) {

      currentTitre_.push_back(initialTitreNormalVariable());

  }

  calcProportions();

}


void kgd::MCMCTITRE::proportion2Titre(const std::vector<double>& proportions) {

  currentTitre_.clear();

  for (auto prop_value : proportions) {

    currentTitre_.push_back(std::log(prop_value));

  }

  proportions_ = proportions;

}


void kgd::MCMCTITRE::updateTitre() {

  std::vector<double> tmpTitre;

  for (size_t k = 0; k < kStrain(); k++) {

    double dt = deltaXnormalVariable();

    currentTitre_[k] += dt;

  }

  calcProportions();

}


void kgd::MCMCTITRE::updateTitreIndex(size_t index) {

  if (index >= kStrain()) {

    ExecEnv::log().critical("MCMCTITRE::updateTitreIndex(); Invalid Index: {} for titre size: {}", index, currentTitre_.size());

  }

  double dt = deltaXnormalVariable();

  currentTitre_[index] += dt;

  calcProportions();

}


double kgd::MCMCTITRE::calcLogPriorTitre() const {

  double log_sum = 0.0;

  for (auto value: currentTitre_) {

    log_sum += std::log(NormalDistribution::pdf(value, mean_log_titre_, sd_log_titre_));

  }

  return log_sum;

}


double kgd::MCMCTITRE::calcPriorTitreIndex(size_t index) const {

  if (index >= kStrain()) {

    ExecEnv::log().critical("MCMCTITRE::calcLogPriorTitreIndex(); Invalid Index: {} for titre size: {}", index, currentTitre_.size());

  }

  return NormalDistribution::pdf(currentTitre_[index], mean_log_titre_, sd_log_titre_);

}


double kgd::MCMCTITRE::deltaXnormalVariable() const {

  return (std_norm_->random(entropy_source_->generator()) * (sd_log_titre_ / update_scale_)) + mean_log_titre_;

}


double kgd::MCMCTITRE::initialTitreNormalVariable() const {

  return (std_norm_->random(entropy_source_->generator()) * sd_log_titre_) + mean_log_titre_;

}


void kgd::MCMCTITRE::calcProportions() {

  proportions_.clear();

  for (auto value : currentTitre_) {

    proportions_.push_back(std::exp(value));

  }

  Utility::normalizeBySum(proportions_);

}

std::string kgd::MCMCTITRE::proportionsText() const {

  std::stringstream ss;
  for (auto value : proportions_) {

    ss << value << ", ";

  }

  return ss.str();

}

std::string kgd::MCMCTITRE::titreText() const {

  std::stringstream ss;
  for (auto value : currentTitre_) {

    ss << value << ", ";

  }

  return ss.str();

}
