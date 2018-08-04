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


// Two implementations of the log beta probability cache.
// The first uses the Armadillo matrix library (SiteProbabilityCache).
// The second just implements using std::vector (stlSiteProbabilityCache).
// This function is a code 'hot spot' and thus an optimization priority.
// Armadillo appears to be about 5% faster for this function.
// This implies that there is no real motivation for converting the rest of the kgd code to Armadillo.

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

  arma::mat llk_surf_mat_;
  arma::rowvec log_beta_gamma_vec_;

  static constexpr size_t A_INDEX = 0;
  static constexpr size_t B_INDEX = 1;


  void armaProbabilityCache(const std::vector<double>& altCount,
                            const std::vector<double>& refCount,
                            double scalingConst,
                            double err,
                            size_t gridSize);

  double armaSiteLogBetaLLK(size_t site, double x) const;

};



class stlSiteProbabilityCache {

public:

  stlSiteProbabilityCache(const std::vector<double>& altCount,
                       const std::vector<double>& refCount,
                       double scalingConst = 100.0,
                       double err = 0.01,
                       size_t gridSize = 99);
  ~stlSiteProbabilityCache() = default;

  double siteLogBetaLLK(size_t site, double x) const { return stlSiteLogBetaLLK(site, x); }

private:

  std::vector<std::vector<double> > llk_surf_;
  std::vector<double> log_beta_gamma_;

  void stlProbabilityCache(const std::vector<double>& altCount,
                           const std::vector<double>& refCount,
                           double scalingConst,
                           double err,
                           size_t gridSize);


  double stlSiteLogBetaLLK(size_t site, double x) const;

};



}   // organization level namespace
}   // project level namespace


#endif //KGD_PROB_CACHE_H
