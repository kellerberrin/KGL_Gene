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
  ~MCMCHAP() =default;

private:


  size_t strainIndex_;
  size_t strainIndex1_;
  size_t strainIndex2_;

  std::shared_ptr<RandomGenerator> mcmcEventRg_;
  double currentLogPriorTitre_;

  /* Structural Subroutines */
  void initializeMcmcChain();

  int sampleMcmcEvent() override;

  void finalizeMcmc() override;

  /* Implementation Subroutines */
  double deltaXnormalVariable() { return stdNorm_->genReal() * SD_LOG_TITRE * 1.0 / PROP_SCALE + MN_LOG_TITRE; }

  double calcLogPriorTitre(std::vector<double> &tmpTitre);

  void writeLastFwdProb(bool useIBD);

  void updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain);

  void initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo);

  void updateProportion();

  std::vector<double> calcTmpTitre();

  double deltaLLKs(std::vector<double> &newLLKs);

  void updateSingleHap();

  void findUpdatingStrainSingle();

  void updatePairHaps();

  void findUpdatingStrainPair();


};


}   // organization level namespace
}   // project level namespace


#endif //KGD_MCMC_HAP_H
