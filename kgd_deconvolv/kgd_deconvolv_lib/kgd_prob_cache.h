//
// Created by kellerberrin on 17/07/18.
//

#ifndef KGL_KGD_PROB_CACHE_H
#define KGL_KGD_PROB_CACHE_H


#include <vector>
#include <iostream>
#include <sstream>
#include <armadillo>
#include "kgd_utility.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace




class SiteProbabilityCache {

public:

  SiteProbabilityCache(const std::vector<double>& altCount,
                       const std::vector<double>& refCount,
                       double scalingConst = 100.0,
                       double err = 0.01,
                       size_t gridSize = 99);
  ~SiteProbabilityCache() = default;

  double siteLogBetaLLK(size_t site, double x) const { return armaSiteLogBetaLLK(site, x); }

private:

  std::vector<std::vector<double> > llk_surf_;
  std::vector<double> log_beta_gamma_;

  arma::mat llk_surf_vec_;
  arma::rowvec log_beta_gamma_vec_;

  void stlProbabilityCache(const std::vector<double>& altCount,
                           const std::vector<double>& refCount,
                           double scalingConst,
                           double err,
                           size_t gridSize);

  void armaProbabilityCache(const std::vector<double>& altCount,
                            const std::vector<double>& refCount,
                            double scalingConst,
                            double err,
                            size_t gridSize);

  double stlSiteLogBetaLLK(size_t site, double x) const;
  double armaSiteLogBetaLLK(size_t site, double x) const;

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_PROB_CACHE_H
