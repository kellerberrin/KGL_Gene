//
// Created by kellerberrin on 30/06/18.
//

#ifndef KGD_CTL_PARAMETER_H
#define KGD_CTL_PARAMETER_H

#include <cstdlib>             // strtol, strtod

namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace

//************************************************************************************
//
// Holds parameters common to the IBD and Haplotype MCMC objects.
//
//************************************************************************************


class MixtureParameterObj {

public:

  explicit MixtureParameterObj();
  MixtureParameterObj(const MixtureParameterObj& copy) = default;
  virtual ~MixtureParameterObj() = default;

  MixtureParameterObj &operator=(const MixtureParameterObj &copy) = default;

  // Get parameters.
  size_t randomSeed() const { return randomSeed_; }

  double parameterG() const { return parameterG_; }
  double getMissCopyProb() const { return missCopyProb_; }
  double constRecombProb() const { return constRecombProb_; }
  double averageCentimorganDistance() const { return averageCentimorganDistance_; }
  double baseCountError() const { return baseCountError_; }

  // Set parameters.
  void setRandomSeed(size_t seed) { randomSeed_ = seed; }

  void setParameterG(double setTo) { parameterG_ = setTo; }
  void setMissCopyProb(double setTo) { missCopyProb_ = setTo; }
  void setConstRecomProb(double setTo) { constRecombProb_ = setTo; }
  void setAvCentiMorgonDistance(double setTo) { averageCentimorganDistance_ = setTo; }
  void setBaseCountError(double setTo) { baseCountError_ = setTo; }


private:


  size_t randomSeed_;  // Random seed the RG

  double baseCountError_; // Probability that a base is miss-counted as reference or alternative.
  double missCopyProb_;  // Probability of miss copy during recombination.
  double averageCentimorganDistance_; // Av recombination distance in centiMorgans.
  double constRecombProb_;  // Recombination scaling factor 1.
  double parameterG_;  // Recombination scaling factor 3.


};


//************************************************************************************
//
// Holds parameters unique to each of the IBD and Haplotype MCMC parameter objects.
//
//************************************************************************************

class MCMCParameterObj : public MixtureParameterObj {

protected:

  MCMCParameterObj() = default;

public:

  MCMCParameterObj(const MCMCParameterObj& copy) = default;
  ~MCMCParameterObj() override = default;

  MCMCParameterObj &operator=(const MCMCParameterObj &copy) = default;

  // Get
  size_t McmcSample() const { return nMcmcSample_; }
  size_t McmcMachineryRate() const { return mcmcMachineryRate_; }
  double McmcBurn() const { return mcmcBurn_; }
  double proposalSigma() const { return parameterSigma_; }
  double proposalMean() const { return hapMean_; }
  double proposalUpdateScaling() const { return scalingFactor_; }
  // Set
  void setMcmcSample(size_t setTo) { nMcmcSample_ = setTo; }
  void setMcmcMachineryRate(size_t setTo) { mcmcMachineryRate_ = setTo; }
  void setMcmcBurn(double setTo) { mcmcBurn_ = setTo; }
  void setProposalSigma(double setTo) { parameterSigma_ = setTo; }
  void setProposalMean(double setTo) { hapMean_ = setTo; }
  void setProposalScaling(double setTo) { scalingFactor_ = setTo; }

private:

  size_t nMcmcSample_;    // Number of MCMC samples (default value 800).
  size_t mcmcMachineryRate_;  // MCMC sample rate (default value 5).
  double mcmcBurn_;  // MCMC burn sample proportion (default value 0.5).

  double parameterSigma_;  // S.D. of the haplotype MCMC proposal
  double scalingFactor_; // haplotype proposal scaling factor.
  double hapMean_; // haplotype proposal mean.

};


//************************************************************************************
//
// The parameter object for the Haplotype MCMC
//
//************************************************************************************


class HapParameterObj : public MCMCParameterObj {

public:

  explicit HapParameterObj();
  HapParameterObj(const HapParameterObj& copy) = default;
  ~HapParameterObj() override = default;

  HapParameterObj &operator=(const HapParameterObj &copy) = default;


private:


};


//************************************************************************************
//
// The parameter object for the IBD MCMC
//
//************************************************************************************

class IBDParameterObj : public MCMCParameterObj {

public:

  explicit IBDParameterObj();
  IBDParameterObj(const IBDParameterObj& copy) = default;
  ~IBDParameterObj() override = default;

  IBDParameterObj &operator=(const IBDParameterObj &copy) = default;

  double betaBinomialConstant() const { return  beta_binomial_constant_; }
  void setBinomialConstant(double setTo) { beta_binomial_constant_ = setTo; }

  size_t cacheGridSize() const { return  cache_grid_size_; }
  void setCacheGridSize(size_t setTo) { cache_grid_size_ = setTo; }

private:


  double beta_binomial_constant_;
  size_t cache_grid_size_;

};






}   // organization level namespace
}   // project level namespace



#endif //KGD_CTL_PARAMETER_H
