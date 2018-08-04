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
#include "kgd_mcmc_titre.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCIBD :public MCMCBASE {

public:

  MCMCIBD(std::shared_ptr<DEploidIO> dEplioidIO,
          std::shared_ptr<McmcSample> mcmcSample);
  ~MCMCIBD() override =default;


private:

  IBDpath ibdPath_;

  size_t theta_accept_;

  /* Structural Subroutines */
  void initializeMcmcChain();

  int sampleMcmcEvent() override;

  void finalizeMcmc() override;

  void recordMcmcMachinery() override;

  /* Implementation Subroutines */
  std::vector<double> computeLlkAtAllSites(const std::vector<double>& proportion, double read_error_prob);

  std::vector<double> averageProportion(const std::vector<std::vector<double> > &proportion);

  void ibdInitializeEssentials(double read_error_prob);

  void ibdSampleMcmcEventStep();

  void ibdUpdateHaplotypesFromPrior();

  void ibdUpdateProportionGivenHap();

  void ibdUpdateThetaGivenHap(double previous_theta, const std::vector<double>& llkAtAllSites);

  void updateProportion();

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_MCMC_IDB_H
