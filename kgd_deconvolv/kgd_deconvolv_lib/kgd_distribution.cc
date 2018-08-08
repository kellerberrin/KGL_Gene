//
// Created by kellerberrin on 5/08/18.
//


#include "kgd_distribution.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp>

namespace bm = boost::math;
namespace kgd = kellerberrin::deconvolv;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Distribution functions use boost special functions.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


double kgd::NormalDistribution::pdf(double x, double mean, double std_dev) {

  static const double inv_sqrt_2pi = 0.3989422804014327;

  double a = (x - mean) / std_dev;

  return (inv_sqrt_2pi / std_dev) * std::exp(-0.5 * a * a);

}

double kgd::BetaDistribution::logInverseBetaFunction(double a, double b) {

  assert(a >= 0);
  assert(b >= 0);

  double log_gamma = bm::lgamma<double>(a + b) - bm::lgamma<double>(a) - bm::lgamma<double>(b);

  return log_gamma;

}


double kgd::BetaDistribution::logPartialPdf(double x, double a, double b) {

  assert(x >= 0 && x <= 1);
  assert(a >= 0);
  assert(b >= 0);

  double ret =  (b - 1) * std::log(1 - x) + (a - 1) * std::log(x);

  return ret;

}


double kgd::BetaDistribution::logPdf(double x, double a, double b) {

  assert(x >= 0 && x <= 1);
  assert(a >= 0);
  assert(b >= 0);

  double ret = logInverseBetaFunction(a, b) + logPartialPdf(x, a, b);

  return ret;

}

double kgd::BetaDistribution::pdf(double x, double a, double b) {

  assert(x >= 0 && x <= 1);
  assert(a >= 0);
  assert(b >= 0);

  double p = bm::tgamma<double>(a + b) / (bm::tgamma<double>(a) * bm::tgamma<double>(b));
  double q = std::pow(1 - x, b - 1) * std::pow(x, a - 1);
  return p * q;

}

double kgd::BetaBinomialDistribution::pdf(size_t n, size_t k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0);
  assert(beta > 0);

  double coeff = bm::binomial_coefficient<double>(n, k);

  double a1 = static_cast<double>(k) + alpha;
  double b1 = static_cast<double>(n-k) + beta;

  double p = bm::beta<double>(a1, b1);
  double q = bm::beta<double>(alpha, beta);

  return coeff * (p / q);

}


double kgd::BetaBinomialDistribution::partialPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0);
  assert(beta > 0);

  double a1 = k + alpha;
  double b1 = n - k + beta;

  double p = bm::beta<double>(a1, b1);
  double q = bm::beta<double>(alpha, beta);

  return (p / q);

}

double kgd::BetaBinomialDistribution::logPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  double r = n - k;
  double a1 = k + alpha;
  double b1 = r + beta;

  double beta_coeff = bm::lgamma<double>(r + 1.0) + bm::lgamma<double>(k + 1.0) - bm::lgamma<double>(n + 2.0);
  double inverse_log_binonimal_coeff = std::log(n + 1.0) + beta_coeff;
  double beta_n = bm::lgamma<double>(a1) + bm::lgamma<double>(b1) - bm::lgamma<double>(a1 + b1);
  double beta_d = bm::lgamma<double>(alpha) + bm::lgamma<double>(beta) - bm::lgamma<double>(alpha + beta);
  double prob = beta_n - beta_d - inverse_log_binonimal_coeff;

  return prob;

}


double kgd::BinomialDistribution::pdf(size_t n, size_t k, double prob_success) {

  assert(k <= n);
  assert(prob_success >= 0);
  assert(prob_success <= 1.0);

  double coeff = bm::binomial_coefficient<double>(n, k);

  double p = std::pow(prob_success, static_cast<double>(k));

  double q = std::pow((1.0 - prob_success), static_cast<double>(n - k));

  return coeff * p * q;

}
