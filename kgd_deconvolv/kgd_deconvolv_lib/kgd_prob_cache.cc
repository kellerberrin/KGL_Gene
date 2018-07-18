//
// Created by kellerberrin on 17/07/18.
//


#include "kgd_prob_cache.h"


namespace kgd = kellerberrin::deconvolv;



kgd::SiteProbabilityCache::SiteProbabilityCache(const std::vector<double>& altCount,
                                                const std::vector<double>& refCount,
                                                double scalingConst,
                                                double err,
                                                size_t gridSize) {

  double pGridSpacing = 1.0 / static_cast<double>(gridSize + 1);

  std::vector<double> pGrid;

  pGrid.push_back(pGridSpacing);

  for (size_t i = 1; i < gridSize; i++) {

    pGrid.push_back(pGrid.back() + pGridSpacing);

  }

  assert(pGrid.size() == gridSize);

  assert(llk_surf_.size() == 0);

  for (size_t i = 0; i < altCount.size(); i++) {

    double alt = altCount[i];
    double ref = refCount[i];

    std::vector<double> ll;

    for (double unadjustedP : pGrid) {

      ll.push_back(Utility::calcLLK(ref, alt, unadjustedP, err, scalingConst));

    }

    double llmax = Utility::max_value(ll);
    std::vector<double> ln;

    for (double lltmp : ll) {

      ln.push_back(std::exp(lltmp - llmax));

    }

    double lnSum = Utility::sumOfVec(ln);
    for (size_t i = 0; i < ln.size(); ++i) {

      ln[i] = ln[i] / lnSum;

    }

    std::vector<double> tmpVec1 = Utility::vecProd(ln, pGrid);
    double mn = Utility::sumOfVec(tmpVec1);
    std::vector<double> pGridSq = Utility::vecProd(pGrid, pGrid);
    std::vector<double> tmpVec2 = Utility::vecProd(ln, pGridSq);
    double vr = Utility::sumOfVec(tmpVec2) - (mn * mn);

    double comm = (mn * ((1.0 - mn) / vr)) - 1.0;
    llk_surf_.push_back(std::vector<double>{mn * comm, (1 - mn) * comm});

  }

  for (auto sitellk : llk_surf_) {

    log_beta_gamma_.push_back(Utility::logBetaGamma(sitellk[0], sitellk[1]));

  }

  assert(llk_surf_.size() == refCount.size());
  assert(log_beta_gamma_.size() == llk_surf_.size());

}


double kgd::SiteProbabilityCache::siteLogBetaLLK(size_t site, double x) const {

  return log_beta_gamma_[site] + Utility::partialLogBetaPdf(x, llk_surf_[site][0], llk_surf_[site][1]);

}
