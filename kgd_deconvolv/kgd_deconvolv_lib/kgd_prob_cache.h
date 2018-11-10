//
// Created by kellerberrin on 17/07/18.
//

#ifndef KGL_KGD_PROB_CACHE_H
#define KGL_KGD_PROB_CACHE_H


#include <vector>
#include <iostream>
#include <sstream>
#include "kgd_utility.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


// Two implementations of the log beta probability cache.
// The first uses the Armadillo matrix library (SiteProbabilityCache).
// The second just implements using std::vector (stlSiteProbabilityCache).
// This function is a code 'hot spot' and thus an optimization priority.
// Armadillo appears to be about 5% faster for this function.
// This implies that there is no real motivation for converting the rest of the kgd code to Armadillo.




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
