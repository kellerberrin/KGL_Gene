//
// Created by kellerberrin on 25/06/18.
//

#ifndef KGD_MCMC_IDB_H
#define KGD_MCMC_IDB_H



#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "randomSample.hpp"   // src/codeCogs/randomSample.hpp
#include "kgd_mcmc_base.h"
#include "kgd_ibdpath.h"
#include "kgd_mcmc_sample.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCIBD :public MCMCBASE {

public:

  MCMCIBD(std::shared_ptr<DEploidIO> dEplioidIO,
          std::shared_ptr<McmcSample> mcmcSample,
          std::shared_ptr<RandomGenerator> randomGenerator);
  ~MCMCIBD() =default;


private:

  IBDpath ibdPath; /* IBD */

  /* Structural Subroutines */
  void initializeMcmcChain();

  int sampleMcmcEvent() override;

  void finalizeMcmc() override;

  /* Implementation Subroutines */

  void initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo);

  std::vector<double> computeLlkAtAllSites(double err = 0.01);

  std::vector<double> averageProportion(const std::vector<std::vector<double> > &proportion);

  void ibdInitializeEssentials();

  void ibdSampleMcmcEventStep();

  void initializePropIBD();

  void ibdUpdateHaplotypesFromPrior();

  void ibdUpdateProportionGivenHap(std::vector<double> &llkAtAllSites);

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_MCMC_IDB_H
