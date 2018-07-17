//
// Created by kellerberrin on 27/06/18.
//

#ifndef KGD_MCMC_TITRE_H
#define KGD_MCMC_TITRE_H

#include <vector>
#include <cstdio>
#include <memory>
#include "randomSample.hpp"

namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCTITRE  {

public:

  MCMCTITRE(size_t strains, double mean_log_titre, double sd_log_titre, double update_scale, size_t random_seed = 0);
  MCMCTITRE(const MCMCTITRE&) = default;
  ~MCMCTITRE() =default;

  MCMCTITRE& operator=(const MCMCTITRE&) = default;

  void randomizeTitre(); // generate an initial titre state.
  void randomizeProportions(); // generate an initial proportion state
  void proportion2Titre(const std::vector<double>& proportions);  // initialize the titre with and proportion vector
  void updateTitre();  // update all the titre elements.
  void updateTitreIndex(size_t index);  // update the indexed element

  double calcLogPriorTitre() const; // calc the sum of the log of the titre probability
  double calcPriorTitreIndex(size_t index) const; // calc the the titre index probability.
  const std::vector<double>& Proportions() const { return proportions_; }
  double hastingsRatio() const { return 1.0; }  // symmetric proposal

  std::string proportionsText() const;
  std::string titreText() const;


private:

  size_t k_strains_;  // number of strains;
  double mean_log_titre_;  // Mean of the norm
  double sd_log_titre_;  // S.D of the norm.
  double update_scale_;  // Update down-scaling

  std::vector<double> currentTitre_;
  std::vector<double> proportions_;

  mutable std::shared_ptr<StandNormalRandomSample> stdNorm_;

  size_t kStrain() const { return k_strains_; }
  double deltaXnormalVariable() const;
  double initialTitreNormalVariable() const;
  void calcProportions();

};



}   // organization level namespace
}   // project level namespace












#endif //KGD_MCMC_TITRE_H
