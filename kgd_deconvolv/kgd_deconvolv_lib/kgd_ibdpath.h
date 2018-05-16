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
namespace deconvolv {          // project level namespace


class IBDpath {
#ifdef UNITTEST
  friend class TestIBDpath;
#endif

public:

  IBDpath() = default;
  ~IBDpath() = default;

  std::vector<std::string> getIBDprobsHeader() const;
  const std::vector<std::vector<double> >& getFwdBwd() const { return fwdbwd_; }
  const std::vector<std::vector<double> >& getLogLikelihoodSurface() const { return llk_surf_; }
  double UpdateHaplotypesFromPrior(size_t strain, size_t loci) const;
  double bestPath(const std::vector<double>& proportion, double err = 0.01) const;

  void init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator);
  void McmcUpdateStep(const std::vector<double>& CurrentProportion);

  void buildPathProbabilityForPainting(const std::vector<double>& proportion); // For painting IBD
  void IBDPathChangeAt(size_t site_idx, double update_divisor) { IBD_path_change_at_[site_idx] /= update_divisor; }

private:

  std::shared_ptr<RandomGenerator> ibd_random_generator_;
  double f_sum_;
  Hprior h_prior_;
  IBDrecombProbs ibd_recomb_probs_;
  std::vector<std::vector<double> > ibd_trans_probs_;
  std::vector<std::vector<double> > fm_;
  std::vector<double> f_sum_state_;
  std::vector<size_t> ibd_configure_path_;

  std::vector<std::vector<double> > bwd_;
  std::vector<std::vector<double> > fwdbwd_;

  size_t kStrain_;
  size_t nLoci_;
  double theta_;

  std::vector<double> current_IBD_path_change_at_;
  std::vector<std::vector<double> > llk_surf_;
  std::vector<int> unique_effectiveK_count_;
  std::vector<double> IBD_path_change_at_;

  void computeAndUpdateTheta();

  void ibdSamplePath(const std::vector<double>& statePrior);

  void computeIbdPathFwdProb(const std::vector<double>& proportion, const std::vector<double>& statePrior);

  std::vector<double> computeStatePrior(const std::vector<double>& effectiveKPrior);

  const Hprior& getHprior() const { return h_prior_; }

  const std::vector<size_t>& idbConfigurePath() const { return ibd_configure_path_; }

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  void setNLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }

  void setTheta(const double setTo) { theta_ = setTo; }

  double theta() const { return theta_; }

  // Methods

  void updateFmAtSiteI(const std::vector<double> &prior, const std::vector<double> &llk);

  void makeIbdTransProbs();

  std::vector<double> computeEffectiveKPrior(double theta);

  void makeLlkSurf(const std::vector<double>& altCount,
                   const std::vector<double>& refCount,
                   double scalingConst = 100.0,
                   double err = 0.01,
                   size_t gridSize = 99);

  void computeUniqueEffectiveKCount();

  std::vector<double> computeLlkOfStatesAtSiteI(const std::vector<double>& proportion, size_t siteI, double err = 0.01);

  std::vector<size_t> findWhichIsSomething(const std::vector<size_t>& tmpOp, size_t something);

  void computeIbdPathBwdProb(const std::vector<double>& proportion, const std::vector<double>& effectiveKPrior, const std::vector<double>& statePrior);

  void combineFwdBwd(const std::vector<std::vector<double>> &reshapedFwd, const std::vector<std::vector<double>> &reshapedBwd);

  std::vector<std::vector<double> > reshapeProbs(const std::vector<std::vector<double> > &probs) const;


};



}   // organization level namespace
}   // project level namespace



#endif //KGD_IBDPATH_H
