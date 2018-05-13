//
// Created by kellerberrin on 13/05/18.
//

#ifndef KGD_IBDPATH_H
#define KGD_IBDPATH_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_mersenne_twister.h"
#include "kgd_deploid_io.h"
#include "kgd_hprior.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace


class IBDpath {
#ifdef UNITTEST
  friend class TestIBDpath;
#endif

  friend class McmcMachinery;
  friend class DEploidIO;

  IBDpath() = default;
  ~IBDpath() = default;

  std::shared_ptr<RandomGenerator> ibdRg_;
  double fSum;
  Hprior hprior;
  IBDrecombProbs ibdRecombProbs;
  std::vector<std::vector<double> > ibdTransProbs;
  std::vector<std::vector<double> > fm;
  std::vector<double> fSumState;
  std::vector<size_t> ibdConfigurePath;

  std::vector<std::vector<double> > bwd;
  std::vector<std::vector<double> > fwdbwd;

  size_t kStrain_;
  size_t nLoci_;
  double theta_;

  std::vector<double> currentIBDpathChangeAt;
  std::vector<std::vector<double> > llkSurf;
  std::vector<int> uniqueEffectiveKCount;
  std::vector<double> IBDpathChangeAt;


  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  void setNLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }

  void setTheta(const double setTo) { theta_ = setTo; }

  double theta() const { return theta_; }


  // Methods
  void computeAndUpdateTheta();

  void updateFmAtSiteI(std::vector<double> &prior,
                       std::vector<double> &llk);

  void ibdSamplePath(std::vector<double> statePrior);

  void makeIbdTransProbs();

  std::vector<double> computeEffectiveKPrior(double theta);

  std::vector<double> computeStatePrior(std::vector<double> effectiveKPrior);

  void makeLlkSurf(std::vector<double> altCount,
                   std::vector<double> refCount,
                   double scalingConst = 100.0,
                   double err = 0.01,
                   size_t gridSize = 99);

  void computeUniqueEffectiveKCount();

  std::vector<double> computeLlkOfStatesAtSiteI(std::vector<double> proportion, size_t siteI, double err = 0.01);

  std::vector<size_t> findWhichIsSomething(std::vector<size_t> tmpOp, size_t something);

  // For painting IBD
  void buildPathProbabilityForPainting(std::vector<double> proportion);

  void computeIbdPathFwdProb(std::vector<double> proportion, std::vector<double> statePrior);

  void computeIbdPathBwdProb(std::vector<double> proportion, std::vector<double> effectiveKPrior, std::vector<double> statePrior);

  void combineFwdBwd(std::vector<std::vector<double>> &reshapedFwd, std::vector<std::vector<double>> &reshapedBwd);

  std::vector<std::vector<double> > reshapeProbs(std::vector<std::vector<double> > &probs);

  double bestPath(std::vector<double> proportion, double err = 0.01);

public:

  std::vector<std::string> getIBDprobsHeader();

  void init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator);

};



}   // organization level namespace
}   // project level namespace



#endif //KGD_IBDPATH_H
