//
// Created by kellerberrin on 25/06/18.
//

#ifndef KGL_KGD_MCMC_BASE_H
#define KGL_KGD_MCMC_BASE_H



#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "kgd_random_generator.h"
#include "randomSample.hpp"   // src/codeCogs/randomSample.hpp
#include "kgd_deploid_io.h"
#include "kgd_panel.h"
#include "kgd_mcmc_sample.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MCMCBASE {

public:

  MCMCBASE(std::shared_ptr<DEploidIO> dEplioidIO,
           std::shared_ptr<McmcSample> mcmcSample,
           std::shared_ptr<RandomGenerator> randomGenerator);

  ~MCMCBASE() = default;

  void runMcmcChain(bool showProgress = true);

protected:


  std::shared_ptr<McmcSample> mcmcSample_;
  std::shared_ptr<DEploidIO> dEploidIO_;
  std::shared_ptr<Panel> panel_;
  size_t seed_;
  std::shared_ptr<RandomGenerator> hapRg_;
  std::shared_ptr<RandomGenerator> propRg_;
  std::shared_ptr<RandomGenerator> initialHapRg_;
  std::shared_ptr<StandNormalRandomSample> stdNorm_;

  std::vector<double> currentProp_;
  std::vector<double> currentLLks_;

  std::vector<std::vector<double> > currentHap_;
  std::vector<double> currentTitre_;

  std::vector<double> currentExpectedWsaf_;
  std::vector<double> cumExpectedWsaf_;

  double MN_LOG_TITRE;
  double SD_LOG_TITRE;
  double PROP_SCALE;

  double burnIn_;
  size_t maxIteration_;
  size_t mcmcThresh_;
  size_t McmcMachineryRate_;

  int eventInt_;

  size_t currentMcmcIteration_;
  bool recordingMcmcBool_;

  virtual int sampleMcmcEvent() = 0;
  virtual void finalizeMcmc() = 0;

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

  void computeDiagnostics();

  void calcMaxIteration(size_t nSample, size_t McmcMachineryRate, double burnIn);

  void recordMcmcMachinery();

  void incrementAccept() { ++acceptUpdate_; }
  int acceptCount() const { return acceptUpdate_; }

private:

  size_t kStrain_;
  size_t nLoci_;

  int acceptUpdate_;

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_MCMC_BASE_H
