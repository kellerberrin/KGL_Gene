//
// Created by kellerberrin on 19/07/18.
//

#ifndef KGL_KGD_IBDPATH_ARMA_H
#define KGL_KGD_IBDPATH_ARMA_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_mersenne_twister.h"
#include "kgd_deploid_io.h"
#include "kgd_hprior.h"
#include "kgd_ibdrecombprobs.h"
#include "kgd_prob_cache.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class IBDPath {
#ifdef UNITTEST
  friend class TestIBDpath;
#endif

public:

  IBDPath() = default;
  ~IBDPath() = default;

  std::vector<std::string> getIBDprobsHeader() const;
  const std::vector<std::vector<double> >& getFwdBwd() const { return fwdbwd_; }
  double siteLogBetaLLK(size_t site, double x) const { return ibd_prob_cache_->siteLogBetaLLK(site, x); }
  double UpdateHaplotypesFromPrior(size_t strain, size_t loci) const;
  double bestPath(const std::vector<double>& proportion, double err = 0.01) const;

  void init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator);
  void McmcUpdateStep(const std::vector<double>& CurrentProportion);

  void buildPathProbabilityForPainting(const std::vector<double>& proportion); // For painting IBD
  void IBDPathChangeAt(size_t site_idx, double update_divisor) { IBD_path_change_at_[site_idx] /= update_divisor; }

private:

  std::shared_ptr<RandomGenerator> ibd_random_generator_;
  std::shared_ptr<const IBDRecombProbs> ibd_recomb_probs_;
  std::shared_ptr<const SiteProbabilityCache> ibd_prob_cache_;
  Hprior h_prior_;

  double f_sum_;
  arma::mat ibd_trans_probs_;
  arma::mat fm_;
  arma::rowvec f_sum_state_;
  arma::Row<size_t> ibd_configure_path_;

  arma::mat bwd_;
  arma::mat fwdbwd_;

  double theta_;
  size_t kStrain_;
  size_t nLoci_;

  arma::rowvec current_IBD_path_change_at_;
  arma::Row<size_t> unique_effectiveK_count_;
  arma::rowvec IBD_path_change_at_;

  void computeAndUpdateTheta();

  void ibdSamplePath(const arma::rowvec& statePrior);

  void computeIbdPathFwdProb(const std::vector<double>& proportion, const arma::rowvec& statePrior);

  void computeStatePrior(arma::rowvec& statePrior, const arma::rowvec& effectiveKPrior);

  const Hprior& getHprior() const { return h_prior_; }

  const std::vector<size_t>& idbConfigurePath() const { return ibd_configure_path_; }

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  void setNLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }

  void setTheta(const double setTo) { theta_ = setTo; }

  double theta() const { return theta_; }

  // Methods

  void updateFmAtSiteI(const arma::rowvec &prior, const arma::rowvec &llk);

  void makeIbdTransProbs();

  void computeEffectiveKPrior(arma::rowvec& effectiveKPrior, double theta);

  void computeUniqueEffectiveKCount();

  void computeLlkOfStatesAtSiteI(arma::rowvec& lk, const std::vector<double>& proportion, size_t siteI, double err = 0.01);

  std::vector<size_t> findWhichIsSomething(const std::vector<size_t>& tmpOp, size_t something);

  void computeIbdPathBwdProb(const std::vector<double>& proportion, const std::vector<double>& effectiveKPrior, const std::vector<double>& statePrior);

  void combineFwdBwd(const std::vector<std::vector<double>> &reshapedFwd, const std::vector<std::vector<double>> &reshapedBwd);

  std::vector<std::vector<double> > reshapeProbs(const std::vector<std::vector<double> > &probs) const;


};



}   // organization level namespace
}   // project level namespace



#endif //KGL_KGD_IBDPATH_ARMA_H
