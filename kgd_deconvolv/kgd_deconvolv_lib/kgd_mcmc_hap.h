//
// Created by kellerberrin on 25/06/18.
//

#ifndef KGD_MCMC_HAP_H
#define KGD_MCMC_HAP_H


#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "kgd_ibdpath.h"
#include "kgd_mcmc_base.h"
#include "kgd_mcmc_sample.h"
#include "kgd_mcmc_titre.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convenience randomization objects.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class HapUpdateType : size_t { UPDATE_PROPORTION = 0, UPDATE_SINGLE_HAPLOTYPE = 1, UPDATE_PAIR_HAPLOTYPE = 2 };
class RandomUpdateType : private RandomInteger {

public:

  RandomUpdateType() : RandomInteger(static_cast<size_t>(HapUpdateType::UPDATE_PROPORTION), static_cast<size_t>(HapUpdateType::UPDATE_PAIR_HAPLOTYPE)) {}
  ~RandomUpdateType() override = default;

  HapUpdateType generateUpdate() const { return static_cast<HapUpdateType>(generate(entropy_source_.generator())); }

private:

  EntropySource entropy_source_;

};


class RandomStrain : private RandomInteger {

public:

  RandomStrain(size_t k_strains) : RandomInteger(0, k_strains - 1) {}
  ~RandomStrain() override = default;

  // Returns a random number 0,..., strain-1
  size_t getRandomStrain() const { return generate(entropy_source_.generator()); }

  // Returns a random pair of different k1 != k2 of strains.
  std::pair<size_t, size_t> getRandomStrainPair() {

    size_t k1 = getRandomStrain();
    size_t k2;
    do {

      k2 = getRandomStrain();

    } while(k1 == k2);

    return {k1, k2};

  }

private:

  EntropySource entropy_source_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MCMC Haplotype object.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class MCMCHAP : public MCMCBASE {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif

public:

  MCMCHAP(std::shared_ptr<DEploidIO> dEplioidIO,
          std::shared_ptr<McmcSample> mcmcSample);
  ~MCMCHAP() override = default;

private:

  RandomUpdateType rand_update_type_;
  RandomStrain rand_strain_;

  std::shared_ptr<Panel> panel_;


  /* Structural Subroutines */
  void initializeMcmcChain();

  int sampleMcmcEvent() override;

  void finalizeMcmc() override;

  void recordMcmcMachinery() override;

  /* Implementation Subroutines */

  void writeLastFwdProb(const std::vector<double>& proportions, bool useIBD);

  void updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain);

  void initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo);

  void updateProportion();

  double deltaLLKs(std::vector<double> &newLLKs);

  void updateSingleHap(const std::vector<double>& proportions);

  void updatePairHaps(const std::vector<double>& proportions);

};


}   // organization level namespace
}   // project level namespace


#endif //KGD_MCMC_HAP_H
