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

  virtual ~MCMCBASE() = default;

protected:


  std::shared_ptr<Panel> panel_;
  std::shared_ptr<RandomGenerator> hapRg_;
  std::shared_ptr<RandomGenerator> propRg_;
  std::shared_ptr<RandomGenerator> initialHapRg_;

  std::vector<double> currentProp_;
  std::vector<double> currentLLks_;

  std::vector<std::vector<double> > currentHap_;
  std::vector<double> currentTitre_;

  std::vector<double> currentExpectedWsaf_;
  std::vector<double> cumExpectedWsaf_;

  double MN_LOG_TITRE;
  double SD_LOG_TITRE;
  double PROP_SCALE;


  void computeDiagnostics() override;
  void recordMcmcMachinery() override;

  void writeLastFwdProb(bool useIBD);

  /* Debug */
  bool doutProp();
  bool doutLLK();

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }
  size_t kStrain() const { return kStrain_; }

  void setNLoci(const size_t setTo) { nLoci_ = setTo; }
  size_t nLoci() const { return nLoci_; }

  void initializeProp();

  void initializeTitre();

  void initializeHap();

  void initializeExpectedWsaf();

  double rBernoulli(double p);

  std::vector<double> titre2prop(std::vector<double> &tmpTitre);

  std::vector<double> calcExpectedWsaf(std::vector<double> &proportion);

  double initialTitreNormalVariable() { return stdNorm_->genReal() * SD_LOG_TITRE + MN_LOG_TITRE; }

private:

  size_t kStrain_;
  size_t nLoci_;

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_MCMC_BASE_H
