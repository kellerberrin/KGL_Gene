//
// Created by kellerberrin on 25/06/18.
//

#ifndef KGD_MCMC_BASE_H
#define KGD_MCMC_BASE_H



#include "kgd_panel.h"
#include "kgd_mcmc_virtual.h"
#include "kgd_ctl_parameter.h"
#include "kgd_mcmc_titre.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCBASE : public MCMCVIRTUAL {

public:

  MCMCBASE(std::shared_ptr<DEploidIO> dEplioidIO,
           std::shared_ptr<McmcSample> mcmcSample,
           const MCMCParameterObj& MCMCParameters);

  ~MCMCBASE() override = default;

protected:

  MCMCParameterObj MCMCParameters_;
  MCMCTITRE titre_proportions_; /* MCMC proportions */

  EntropySource entropy_source_;  // Entropy for the random generators.
  UniformUnitDistribution random_unit_;  // Uniform [0, 1] random generator;

  std::vector<std::vector<double> > currentHap_;
  std::vector<double> currentExpectedWsaf_;
  std::vector<double> cumExpectedWsaf_;
  std::vector<double> currentLLks_;

  void computeDiagnostics() override;

  void writeLastFwdProb(bool useIBD);

  /* Debug */

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }
  size_t kStrain() const { return kStrain_; }

  void setNLoci(const size_t setTo) { nLoci_ = setTo; }
  size_t nLoci() const { return nLoci_; }

  void initializeHap();

  void initializeExpectedWsaf(const std::vector<double>& proportion);

  double rBernoulli(double p);

  std::vector<double> calcExpectedWsaf(const std::vector<double>& proportion);

private:

  size_t kStrain_;
  size_t nLoci_;

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_MCMC_BASE_H
