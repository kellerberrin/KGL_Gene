//
// Created by kellerberrin on 25/06/18.
//

#ifndef KGD_MCMC_BASE_H
#define KGD_MCMC_BASE_H



#include "kgd_panel.h"
#include "kgd_mcmc_virtual.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCBASE : public MCMCVIRTUAL {

public:

  MCMCBASE(std::shared_ptr<DEploidIO> dEplioidIO,
           std::shared_ptr<McmcSample> mcmcSample,
           std::shared_ptr<RandomGenerator> randomGenerator);

  ~MCMCBASE() override = default;

protected:


  std::shared_ptr<RandomGenerator> hapRg_;
  std::shared_ptr<RandomGenerator> propRg_;
  std::shared_ptr<RandomGenerator> initialHapRg_;

  std::vector<double> currentLLks_;

  std::vector<std::vector<double> > currentHap_;

  std::vector<double> currentExpectedWsaf_;
  std::vector<double> cumExpectedWsaf_;

  void computeDiagnostics() override;

  void writeLastFwdProb(bool useIBD);

  /* Debug */
  bool doutLLK();

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
