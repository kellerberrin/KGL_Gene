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


class MCMCHAP : public MCMCBASE {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif

public:

  MCMCHAP(std::shared_ptr<DEploidIO> dEplioidIO,
          std::shared_ptr<McmcSample> mcmcSample,
          std::shared_ptr<RandomGenerator> randomGenerator);
  ~MCMCHAP() override = default;

private:


  size_t strainIndex_;
  size_t strainIndex1_;
  size_t strainIndex2_;

  std::shared_ptr<RandomGenerator> mcmcEventRg_;
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

  void findUpdatingStrainSingle();

  void updatePairHaps(const std::vector<double>& proportions);

  void findUpdatingStrainPair();



};


}   // organization level namespace
}   // project level namespace


#endif //KGD_MCMC_HAP_H
