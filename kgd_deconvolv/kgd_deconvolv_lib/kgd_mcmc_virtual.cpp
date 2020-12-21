//
// Created by kellerberrin on 26/06/18.
//


#include <cstdio>
#include <cmath>       // ceil
#include "kgd_deconvolv_app.h"
#include "kgd_utility.h"
#include "kgd_mcmc_virtual.h"

namespace kgd = kellerberrin::deconvolv;


kgd::MCMCVIRTUAL::MCMCVIRTUAL(std::shared_ptr<DEploidIO> dEploidIO,
                              std::shared_ptr<McmcSample> mcmcSample) {

  dEploidIO_ = dEploidIO;
  mcmcSample_ = mcmcSample;
  acceptUpdate_ = 0;
  record_count_ = 0;

}


void kgd::MCMCVIRTUAL::runMcmcChain(bool showProgress) {

  for (currentMcmcIteration_ = 1; currentMcmcIteration_ <= maxIteration_; ++currentMcmcIteration_) {

    eventInt_ = sampleMcmcEvent();

    if (currentMcmcIteration_ % 100 == 0 && showProgress) {

      ExecEnv::log().info("MCMC iteration: {}/{}, completed: {}%", currentMcmcIteration_, maxIteration_, int(currentMcmcIteration_ * 100 / maxIteration_));

    }

    recordingMcmcBool_ = (currentMcmcIteration_ > mcmcThresh_ && currentMcmcIteration_ % McmcMachineryRate_ == 0);

    if (recordingMcmcBool_) {

      ++record_count_;
      recordMcmcMachinery();

    }

  }

  finalizeMcmc();

  computeDiagnostics();

  ExecEnv::log().info("#### MCMC RUN finished ####");

}


void kgd::MCMCVIRTUAL::calcMaxIteration(size_t nSample, size_t McmcMachineryRate, double burnIn) {

  burnIn_ = burnIn;  // Burn in proportion
  McmcMachineryRate_ = McmcMachineryRate;  // Thinning rate

  // MCMC iterations including burnin.
  double fmax_iter = static_cast<double>(nSample) * static_cast<double>(McmcMachineryRate) / (1.0 - burnIn_);
  maxIteration_ = static_cast<size_t>(std::ceil(fmax_iter));

  // MCMC iterations after burnin
  double fmcmc_thresh = static_cast<double>(nSample) * static_cast<double>(McmcMachineryRate) * burnIn_ / (1.0 - burnIn_);
  mcmcThresh_ = static_cast<size_t>(std::ceil(fmcmc_thresh));

}



