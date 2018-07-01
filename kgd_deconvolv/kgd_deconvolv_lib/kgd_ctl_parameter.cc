//
// Created by kellerberrin on 30/06/18.
//

#include "kgd_ctl_parameter.h"


namespace kgd = kellerberrin::deconvolv;


kgd::MixtureParameterObj::MixtureParameterObj() {


  nMcmcSample_ = 800;
  mcmcBurn_ = 0.5;
  mcmcMachineryRate_ = 5;

  missCopyProb_ = 0.01;
  constRecombProb_ = 1.0;
  averageCentimorganDistance_ = 15000.0;

  setScalingFactor(100.0);
  setParameterG(20.0);
  setParameterSigma(5.0);
  setIBDSigma(20.0);
  set_seed(0);
  setKstrain(5);

}
