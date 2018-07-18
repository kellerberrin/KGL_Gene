//
// Created by kellerberrin on 17/07/18.
//

#ifndef KGL_KGD_PROB_CACHE_H
#define KGL_KGD_PROB_CACHE_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
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

  double sitelogBetaLLK(size_t site, double x) const;
  const std::vector<std::vector<double> >& getLogLikelihoodSurface() const { return llk_surf_; }

private:

  std::vector<std::vector<double> > llk_surf_;
  std::vector<double> logBetaGamma_;

};



}   // organization level namespace
}   // project level namespace


#endif //KGL_KGD_PROB_CACHE_H
