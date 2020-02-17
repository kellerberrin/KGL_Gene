//
// Created by kellerberrin on 17/07/18.
//

#include "kgd_deconvolv_app.h"
#include "kel_distribution.h"
#include "kgd_prob_cache.h"


namespace kgd = kellerberrin::deconvolv;



//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************
//*******************************************************************************************


kgd::stlSiteProbabilityCache::stlSiteProbabilityCache(const std::vector<double>& altCount,
                                                const std::vector<double>& refCount,
                                                double scalingConst,
                                                double err,
                                                size_t gridSize) {

  stlProbabilityCache(altCount, refCount, scalingConst, err, gridSize);

}


void kgd::stlSiteProbabilityCache::stlProbabilityCache(const std::vector<double>& altCount,
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

  for (size_t site = 0; site < altCount.size(); ++site) {

    double alt = altCount[site];
    double ref = refCount[site];

    std::vector<double> ll;

    for (double unadjustedP : pGrid) {

      double llk = Utility::calcLLK(ref, alt, unadjustedP, err, scalingConst);

      ll.push_back(llk);

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

    log_beta_gamma_.push_back(BetaDistribution::logInverseBetaFunction(sitellk[0], sitellk[1]));

  }

  assert(llk_surf_.size() == refCount.size());
  assert(log_beta_gamma_.size() == llk_surf_.size());

}

double kgd::stlSiteProbabilityCache::stlSiteLogBetaLLK(size_t site, double x) const {

  return log_beta_gamma_[site] + BetaDistribution::logPartialPdf(x, llk_surf_[site][0], llk_surf_[site][1]);

}
