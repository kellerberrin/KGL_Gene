//
// Created by kellerberrin on 26/06/18.
//

#ifndef KGD_MCMC_VIRT_H
#define KGD_MCMC_VIRT_H


#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "kgd_random.h"
#include "randomSample.hpp"   // src/codeCogs/randomSample.hpp
#include "kgd_deploid_io.h"
#include "kgd_mcmc_sample.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCVIRTUAL {

public:

  MCMCVIRTUAL(std::shared_ptr<DEploidIO> dEplioidIO,
              std::shared_ptr<McmcSample> mcmcSample);

  virtual ~MCMCVIRTUAL() = default;

  void runMcmcChain(bool showProgress = true);

protected:

  virtual int sampleMcmcEvent() = 0;
  virtual void finalizeMcmc() = 0;
  virtual void recordMcmcMachinery() = 0;
  virtual void computeDiagnostics() = 0;

  std::shared_ptr<McmcSample> mcmcSample_;
  std::shared_ptr<DEploidIO> dEploidIO_;

  void calcMaxIteration(size_t nSample, size_t McmcMachineryRate, double burnIn);

  void incrementAccept() { ++acceptUpdate_; }
  size_t acceptCount() const { return acceptUpdate_; }

  size_t total_MCMC_iterations() const { return maxIteration_; }
  size_t burnin_MCMC_iterations() const { return mcmcThresh_; }
  size_t current_MCMC_iteration() const { return currentMcmcIteration_; }

  int eventType() const { return eventInt_; }

private:

  size_t acceptUpdate_;
  int eventInt_;

  size_t currentMcmcIteration_;
  bool recordingMcmcBool_;

  double burnIn_;
  size_t maxIteration_;
  size_t mcmcThresh_;
  size_t McmcMachineryRate_;

};


}   // organization level namespace
}   // project level namespace



#endif //KGD_MCMC_VIRT_H
