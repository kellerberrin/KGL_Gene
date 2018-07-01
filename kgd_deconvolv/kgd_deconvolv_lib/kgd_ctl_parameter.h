//
// Created by kellerberrin on 30/06/18.
//

#ifndef KGD_CTL_PARAMETER_H
#define KGD_CTL_PARAMETER_H

#include <cstdlib>             // strtol, strtod

namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


class MixtureParameterObj {

public:

  explicit MixtureParameterObj();
  ~MixtureParameterObj() = default;

  MixtureParameterObj &operator=(const MixtureParameterObj &copy) = default;

  // Get parameters.
  size_t kStrain() const { return kStrain_; }
  size_t getMcmcSample() const { return nMcmcSample_; }
  size_t getMcmcMachineryRate() const { return mcmcMachineryRate_; }
  double getMcmcBurn() const { return mcmcBurn_; }
  double scalingFactor() const { return scalingFactor_; }
  double getMissCopyProb() const { return missCopyProb_; }
  double ibdSigma() const { return ibdSigma_; }
  double parameterSigma() const { return parameterSigma_; }
  double averageCentimorganDistance() const { return averageCentimorganDistance_; }
  double parameterG() const { return parameterG_; }
  double constRecombProb() const { return constRecombProb_; }
  size_t randomSeed() const { return randomSeed_; }

  // Set parameters.
  void set_seed(const size_t seed) { randomSeed_ = seed; }
  void setParameterG(const double setTo) { parameterG_ = setTo; }
  void setParameterSigma(const double setTo) { parameterSigma_ = setTo; }
  void setIBDSigma(const double setTo) { ibdSigma_ = setTo; }
  void setKstrain(const size_t setTo) { kStrain_ = setTo; }
  void setScalingFactor(const double setTo) { scalingFactor_ = setTo; }


private:


  size_t nMcmcSample_;    // Number of MCMC samples (default value 800).
  size_t mcmcMachineryRate_;  // MCMC sample rate (default value 5).
  double mcmcBurn_;  // MCMC burn rate (default value 0.5).

  size_t kStrain_;
  size_t randomSeed_;

  // Parameters
  double missCopyProb_;
  double averageCentimorganDistance_;// = 15000.0,
  double constRecombProb_;
  double scalingFactor_; // 100.0

  double parameterG_;
  double parameterSigma_;
  double ibdSigma_;


};



}   // organization level namespace
}   // project level namespace



#endif //KGD_CTL_PARAMETER_H
