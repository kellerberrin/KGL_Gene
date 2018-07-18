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

  armaProbabilityCache(altCount, refCount, scalingConst, err, gridSize);

}


void kgd::SiteProbabilityCache::stlProbabilityCache(const std::vector<double>& altCount,
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


void kgd::SiteProbabilityCache::armaProbabilityCache(const std::vector<double>& altCount,
                                                     const std::vector<double>& refCount,
                                                     double scalingConst,
                                                     double err,
                                                     size_t gridSize) {

  arma::rowvec ref_count(refCount);
  arma::rowvec alt_count(altCount);

  arma::rowvec grid(gridSize);

  double pGridSpacing = 1.0 / static_cast<double>(gridSize + 1);

  for (size_t i = 0; i < gridSize; ++i) {

    grid(i) = static_cast<double>(i + 1) * pGridSpacing;

  }

  llk_surf_vec_.resize(alt_count.n_cols, 2);
  log_beta_gamma_vec_.resize(alt_count.n_cols);

  for (size_t site = 0; site < alt_count.n_cols; ++site) {

    arma::rowvec ll(gridSize);

    for (size_t i = 0; i < gridSize; ++i) {

      ll(i) = Utility::calcLLK(ref_count(site), alt_count(site), grid(i), err, scalingConst);

    }

    double llmax = ll.max();

    arma::rowvec ln(gridSize);

    for (size_t i = 0; i < gridSize; ++i) {

      ln(i) = std::exp(ll(i) - llmax);

    }

    double lnSum = arma::accu(ln);
    for (size_t i = 0; i < ln.size(); ++i) {

      ln(i) = ln(i) / lnSum;

    }

    arma::rowvec mn_vec(ln % grid);

    double mn = arma::accu(mn_vec);

    double vr = arma::accu(mn_vec % grid) - (mn * mn);

    double comm = (mn * ((1.0 - mn) / vr)) - 1.0;

    double a = mn * comm;
    double b = (1 - mn) * comm;

    llk_surf_vec_(site,0) = a;
    llk_surf_vec_(site,1) = b;

    log_beta_gamma_vec_(site) = Utility::logBetaGamma(a, b);

  }

}


double kgd::SiteProbabilityCache::stlSiteLogBetaLLK(size_t site, double x) const {

  return log_beta_gamma_[site] + Utility::partialLogBetaPdf(x, llk_surf_[site][0], llk_surf_[site][1]);

}

double kgd::SiteProbabilityCache::armaSiteLogBetaLLK(size_t site, double x) const {

  return log_beta_gamma_vec_(site) + Utility::partialLogBetaPdf(x, llk_surf_vec_(site, 0), llk_surf_vec_(site, 1));

}

