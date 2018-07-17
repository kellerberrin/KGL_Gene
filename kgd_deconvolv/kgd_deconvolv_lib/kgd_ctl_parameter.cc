//
// Created by kellerberrin on 30/06/18.
//

#include "kgd_ctl_parameter.h"


namespace kgd = kellerberrin::deconvolv;




kgd::MixtureParameterObj::MixtureParameterObj() {

  setRandomSeed(0);

  setMissCopyProb(0.01);
  setConstRecomProb(1.0);
  setAvCentiMorgonDistance(15000.0);
  setParameterG(20.0);
  setBaseCountError(0.01);

}


kgd::HapParameterObj::HapParameterObj() {

  setMcmcSample(800);
  setMcmcMachineryRate(5);
  setMcmcBurn(0.5);

  setProposalSigma(5.0);
  setProposalMean(0.0);
  setProposalScaling(100.0);

}


kgd::IBDParameterObj::IBDParameterObj() {

  setMcmcSample(100);
  setMcmcMachineryRate(10);
  setMcmcBurn(0.5);

  setProposalSigma(10.0);  // Originally 20.0
  setProposalMean(0.0);
  setProposalScaling(40.0);

}
