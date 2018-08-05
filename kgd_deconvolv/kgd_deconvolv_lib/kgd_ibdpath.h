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
#include "kgd_deploid_io.h"
#include "kgd_hprior.h"
#include "kgd_ibdrecombprobs.h"
#include "kgd_prob_cache.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class IBDpath {
#ifdef UNITTEST
  friend class TestIBDpath;
#endif

public:

  IBDpath() = default;
  ~IBDpath() = default;

  void init(DEploidIO &dEploidIO);
  void McmcUpdateStep(const std::vector<double>& CurrentProportion);

  void buildPathProbabilityForPainting(const std::vector<double>& proportion); // For painting IBD
  void IBDPathChangeAt(size_t site_idx, double update_divisor) { IBD_path_change_at_[site_idx] /= update_divisor; }

  std::vector<std::string> getIBDprobsHeader() const;
  double UpdateHaplotypesFromPrior(size_t strain, size_t loci) const;
  double bestPath(const std::vector<double>& proportion, double err = 0.01) const;
  const std::vector<std::vector<double> >& getFwdBwd() const { return fwdbwd_; }
  double siteLogBetaLLK(size_t site, double x) const { return ibd_prob_cache_->siteLogBetaLLK(site, x); }

  void setTheta(const double setTo) { theta_ = setTo; }
  double theta() const { return theta_; }

private:

  std::shared_ptr<const IBDRecombProbs> ibd_recomb_probs_;
  std::shared_ptr<const SiteProbabilityCache> ibd_prob_cache_;
  Hprior h_prior_;

  double f_sum_;
  std::vector<std::vector<double> > fm_;
  std::vector<double> f_sum_state_;

  std::vector<std::vector<double> > ibd_trans_probs_;
  std::vector<size_t> ibd_configure_path_;

  std::vector<std::vector<double> > bwd_;
  std::vector<std::vector<double> > fwdbwd_;

  size_t kStrain_;
  size_t nLoci_;
  double theta_;
  double recomb_miss_copy_prob_;

  std::vector<double> current_IBD_path_change_at_;
  std::vector<size_t> unique_effectiveK_count_;
  std::vector<double> IBD_path_change_at_;

  void computeAndUpdateTheta();

  void ibdSamplePath(const std::vector<double>& statePrior);

  void computeIbdPathFwdProb(const std::vector<double>& proportion, const std::vector<double>& state_prior);

  std::vector<double> computeStatePrior(const std::vector<double>& effectiveKPrior);

  const Hprior& getHprior() const { return h_prior_; }

  const std::vector<size_t>& idbConfigurePath() const { return ibd_configure_path_; }

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  void setNLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }


  // Methods



  void updateFmAtSiteI(const std::vector<double> &prior, const std::vector<double> &llk);

  void makeIbdTransProbs();

  std::vector<double> computeEffectiveKPrior(double theta);

  void computeUniqueEffectiveKCount();

  std::vector<double> computeLlkOfStatesAtSiteI(const std::vector<double>& proportion, size_t site_i, double miss_copy_error);

  std::vector<size_t> findAllIndex(const std::vector<size_t> &index_array, size_t index);

  void computeIbdPathBwdProb(const std::vector<double>& proportion, const std::vector<double>& effectiveKPrior, const std::vector<double>& statePrior);

  void combineFwdBwd(const std::vector<std::vector<double>> &reshapedFwd, const std::vector<std::vector<double>> &reshapedBwd);

  std::vector<std::vector<double> > reshapeProbs(const std::vector<std::vector<double> > &probs) const;


};



}   // organization level namespace
}   // project level namespace



#endif //KGD_IBDPATH_H
