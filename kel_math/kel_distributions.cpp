// Kellerberrin 2026.

#include "kel_exec_env.h"
#include "kel_continuous_distributions.h"
#include "kel_discrete_distributions.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <numeric>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

namespace bm = boost::math;
namespace kel = kellerberrin;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Distribution functions use boost special functions.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Normal distribution
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::NormalDistribution::pdf(double x, double mean, double std_dev) {

  static constexpr double INV_SQRT_2PI = 0.3989422804014327;

  if (std_dev <= 0.0) {

    ExecEnv::log().error("NormalDistribution::pdf; std_dev must be > 0.0, got: {}", std_dev);
    return 0.0;

  }

  const double a = (x - mean) / std_dev;
  return (INV_SQRT_2PI / std_dev) * std::exp(-0.5 * a * a);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gamma distribution
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::GammaDistribution::pdf(double x) const {

  try {

    return bm::pdf(bm::gamma_distribution<>(shape_, scale_), x);

  } catch (const std::exception& e) {

    ExecEnv::log().error("GammaDistribution::pdf; boost error: {}", e.what());
    return 0.0;

  }

}

double kel::GammaDistribution::cdf(double x) const {

  try {

    return bm::cdf(bm::gamma_distribution<>(shape_, scale_), x);

  } catch (const std::exception& e) {

    ExecEnv::log().error("GammaDistribution::cdf; boost error: {}", e.what());
    return 0.0;

  }

}

double kel::GammaDistribution::quantile(double p) const {

  try {

    return bm::quantile(bm::gamma_distribution<>(shape_, scale_), p);

  } catch (const std::exception& e) {

    ExecEnv::log().error("GammaDistribution::quantile; boost error: {}", e.what());
    return 0.0;

  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta distribution
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::BetaDistribution::logInverseBetaFunction(double a, double b) {

  assert(a > 0.0);
  assert(b > 0.0);

  try {

    return bm::lgamma<double>(a + b) - bm::lgamma<double>(a) - bm::lgamma<double>(b);

  } catch (const std::exception& e) {

    ExecEnv::log().error("BetaDistribution::logInverseBetaFunction; boost error: {}", e.what());
    return 0.0;

  }

}


double kel::BetaDistribution::logPartialPdf(double x, double a, double b) {

  assert(x >= 0.0 && x <= 1.0);
  assert(a >= 0.0);
  assert(b >= 0.0);

  return (b - 1.0) * std::log(1.0 - x) + (a - 1.0) * std::log(x);

}


double kel::BetaDistribution::logPdf(double x, double a, double b) {

  assert(x >= 0.0 && x <= 1.0);
  assert(a > 0.0);
  assert(b > 0.0);

  return logInverseBetaFunction(a, b) + logPartialPdf(x, a, b);

}

double kel::BetaDistribution::pdf(double x, double a, double b) {

  assert(x >= 0.0 && x <= 1.0);
  assert(a > 0.0);
  assert(b > 0.0);

  try {

    const double p = bm::tgamma<double>(a + b) / (bm::tgamma<double>(a) * bm::tgamma<double>(b));
    const double q = std::pow(1.0 - x, b - 1.0) * std::pow(x, a - 1.0);
    return p * q;

  } catch (const std::exception& e) {

    ExecEnv::log().error("BetaDistribution::pdf; boost error: {}", e.what());
    return 0.0;

  }

}


double kel::BetaDistribution::mean(double a, double b) {

  assert(a > 0.0);
  assert(b > 0.0);

  return a / (a + b);

}


double kel::BetaDistribution::var(double a, double b) {

  assert(a > 0.0);
  assert(b > 0.0);

  return (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

}


double kel::BetaDistribution::mode(double a, double b) {

  assert(a > 1.0);
  assert(b > 1.0);

  return (a - 1.0) / (a + b - 2.0);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta Binomial distribution
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::BetaBinomialDistribution::pdf(std::size_t n, std::size_t k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  try {

    const double coeff = bm::binomial_coefficient<double>(static_cast<double>(n), static_cast<double>(k));

    const double a1 = static_cast<double>(k) + alpha;
    const double b1 = static_cast<double>(n - k) + beta;

    const double p = bm::beta<double>(a1, b1);
    const double q = bm::beta<double>(alpha, beta);

    return coeff * (p / q);

  } catch (const std::exception& e) {

    ExecEnv::log().error("BetaBinomialDistribution::pdf; boost error: {}", e.what());
    return 0.0;

  }

}


double kel::BetaBinomialDistribution::partialPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  try {

    const double a1 = k + alpha;
    const double b1 = n - k + beta;

    const double p = bm::beta<double>(a1, b1);
    const double q = bm::beta<double>(alpha, beta);

    return p / q;

  } catch (const std::exception& e) {

    ExecEnv::log().error("BetaBinomialDistribution::partialPdf; boost error: {}", e.what());
    return 0.0;

  }

}


double kel::BetaBinomialDistribution::logPartialPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  try {

    const double r = n - k;
    const double a1 = k + alpha;
    const double b1 = r + beta;

    const double beta_n = bm::lgamma<double>(a1) + bm::lgamma<double>(b1) - bm::lgamma<double>(a1 + b1);
    const double beta_d = bm::lgamma<double>(alpha) + bm::lgamma<double>(beta) - bm::lgamma<double>(alpha + beta);

    return beta_n - beta_d;

  } catch (const std::exception& e) {

    ExecEnv::log().error("BetaBinomialDistribution::logPartialPdf; boost error: {}", e.what());
    return 0.0;

  }

}


double kel::BetaBinomialDistribution::logPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  try {

    const double r = n - k;

    // Log of the binomial coefficient C(n,k) = n! / (k! * (n-k)!).
    const double log_binomial_coeff = bm::lgamma<double>(n + 1.0)
                                      - bm::lgamma<double>(k + 1.0)
                                      - bm::lgamma<double>(r + 1.0);

    return logPartialPdf(n, k, alpha, beta) + log_binomial_coeff;

  } catch (const std::exception& e) {

    ExecEnv::log().error("BetaBinomialDistribution::logPdf; boost error: {}", e.what());
    return 0.0;

  }

}


// .first is alpha, .second is beta. The raw moments are calculated from the observations and used to calculate alpha and beta.
std::pair<double, double> kel::BetaBinomialDistribution::methodOfMoments(const std::vector<std::size_t>& observations, std::size_t n_trials) {

  if (observations.empty()) {

    ExecEnv::log().warn("BetaBinomialDistribution::methodOfMoments; empty observations vector, returning (0, 0)");
    return {0.0, 0.0};

  }

  // calculate the first moment
  const std::size_t obs_sum = std::accumulate(observations.begin(), observations.end(), std::size_t{0});
  const double m1 = static_cast<double>(obs_sum) / static_cast<double>(observations.size());

  // calculate the 2nd moment
  const std::size_t sqr_sum = std::transform_reduce(observations.begin(), observations.end(),
                                                    std::size_t{0}, std::plus{},
                                                    [](std::size_t obs) { return obs * obs; });
  const double m2 = static_cast<double>(sqr_sum) / static_cast<double>(observations.size());

  const double n = static_cast<double>(n_trials);

  const double a_numer = (n * m1) - m2;
  const double ab_denom = n * ((m2 / m1) - m1 - 1.0) + m1;

  if (m1 == 0.0 || ab_denom == 0.0) {

    ExecEnv::log().warn("BetaBinomialDistribution::methodOfMoments; degenerate moments, returning (0, 0)");
    return {0.0, 0.0};

  }

  const double a = a_numer / ab_denom;

  const double b_numer = (n - m1) * (n - (m2 / m1));
  const double b = b_numer / ab_denom;

  return {a, b};

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Binomial distribution
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::BinomialDistribution::pdf(std::size_t n, std::size_t k, double prob_success) {

  assert(k <= n);
  assert(prob_success >= 0.0 && prob_success <= 1.0);

  try {

    const double coeff = bm::binomial_coefficient<double>(static_cast<double>(n), static_cast<double>(k));

    const double p = std::pow(prob_success, static_cast<double>(k));
    const double q = std::pow((1.0 - prob_success), static_cast<double>(n - k));

    return coeff * p * q;

  } catch (const std::exception& e) {

    ExecEnv::log().error("BinomialDistribution::pdf; boost error: {}", e.what());
    return 0.0;

  }

}


double kel::BinomialDistribution::cdf(std::size_t n, double k, double prob_success) {

  assert(k <= n);
  assert(prob_success >= 0.0 && prob_success <= 1.0);

  if (k < 0.0) {

    ExecEnv::log().error("BinomialDistribution::cdf; k:{} must be >= 0.0", k);
    return 0.0;

  }

  const std::size_t integer_k = static_cast<std::size_t>(std::floor(k));
  double sum_pdf = 0.0;

  for (std::size_t index = 0; index <= integer_k; ++index) {

    sum_pdf += pdf(n, index, prob_success);

  }

  return sum_pdf;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hypergeometric distribution
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

kel::HypergeometricDistribution::HypergeometricDistribution(std::size_t pop_successes_K,
                                                            std::size_t sample_size_n,
                                                            std::size_t population_N) {

  if (pop_successes_K > population_N) {

    ExecEnv::log().warn("HypergeometricDistribution::HypergeometricDistribution; Population Successes K:{} exceeds population size :{}",
                         pop_successes_K, population_N);
    pop_successes_K = population_N;
  }

  if (sample_size_n > population_N) {

    ExecEnv::log().warn("HypergeometricDistribution::HypergeometricDistribution; Sample size n:{} exceeds population size :{}",
                         sample_size_n, population_N);
    sample_size_n = population_N;
  }

  pop_successes_K_ = pop_successes_K;
  sample_size_n_ = sample_size_n;
  population_N_ = population_N;

}

std::size_t kel::HypergeometricDistribution::clampSuccesses_k(std::size_t successes_k) const {

  if (successes_k > upperSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution::clampSuccesses_k; r_successes k:{} exceeds upper limit :{}",
                        successes_k, upperSuccesses_k());
    successes_k = upperSuccesses_k();

  }

  if (successes_k < lowerSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution::clampSuccesses_k; r_successes k:{} below lower limit :{}",
                        successes_k, lowerSuccesses_k());
    successes_k = lowerSuccesses_k();

  }

  return successes_k;

}

double kel::HypergeometricDistribution::pdf(std::size_t successes_k) const {

  successes_k = clampSuccesses_k(successes_k);

  try {

    bm::hypergeometric_distribution hypergeometric(pop_successes_K_, sample_size_n_, population_N_);
    return bm::pdf(hypergeometric, successes_k);

  } catch (const std::exception& e) {

    ExecEnv::log().error("HypergeometricDistribution::pdf; boost error: {}", e.what());
    return 0.0;

  }

}

double kel::HypergeometricDistribution::cdf(std::size_t successes_k) const {

  successes_k = clampSuccesses_k(successes_k);

  try {

    bm::hypergeometric_distribution hypergeometric(pop_successes_K_, sample_size_n_, population_N_);
    return bm::cdf(hypergeometric, successes_k);

  } catch (const std::exception& e) {

    ExecEnv::log().error("HypergeometricDistribution::cdf; boost error: {}", e.what());
    return 0.0;

  }

}

std::size_t kel::HypergeometricDistribution::quantile(double p) const {

  if (p < 0.0 || p > 1.0) {

    ExecEnv::log().error("HypergeometricDistribution::quantile; probability p:{} must be in [0.0, 1.0]", p);
    return 0;

  }

  try {

    bm::hypergeometric_distribution hypergeometric(pop_successes_K_, sample_size_n_, population_N_);
    return static_cast<std::size_t>(std::llround(bm::quantile(hypergeometric, p)));

  } catch (const std::exception& e) {

    ExecEnv::log().error("HypergeometricDistribution::quantile; boost error: {}", e.what());
    return 0;

  }

}


std::size_t kel::HypergeometricDistribution::lowerSuccesses_k() const {

  const std::int64_t lower_success = static_cast<std::int64_t>(sample_size_n_ + pop_successes_K_) - static_cast<std::int64_t>(population_N_);
  return static_cast<std::size_t>(std::max<std::int64_t>(0, lower_success));

}


double kel::HypergeometricDistribution::upperSingleTailTest(std::size_t test_value_k) const {

  if (test_value_k > 0) {

    --test_value_k;

  }

  return 1.0 - cdf(test_value_k);

}


double kel::HypergeometricDistribution::lowerSingleTailTest(std::size_t test_value_k) const {

  return cdf(test_value_k);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Poisson distribution. Uses boost for implementation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::Poisson::pdf(std::size_t count) const {

  try {

    return bm::pdf(bm::poisson_distribution<>(lambda_), static_cast<double>(count));

  } catch (const std::exception& e) {

    ExecEnv::log().error("Poisson::pdf; boost error: {}", e.what());
    return 0.0;

  }

}

double kel::Poisson::cdf(std::size_t count) const {

  try {

    return bm::cdf(bm::poisson_distribution<>(lambda_), static_cast<double>(count));

  } catch (const std::exception& e) {

    ExecEnv::log().error("Poisson::cdf; boost error: {}", e.what());
    return 0.0;

  }

}

std::size_t kel::Poisson::quantile(double p) const {

  try {

    return static_cast<std::size_t>(bm::quantile(bm::poisson_distribution<>(lambda_), p));

  } catch (const std::exception& e) {

    ExecEnv::log().error("Poisson::quantile; boost error: {}", e.what());
    return 0;

  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Negative Binomial distribution. Uses boost for implementation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::NegativeBinomial::pdf(std::size_t count) const {

  try {

    return bm::pdf(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_), static_cast<double>(count));

  } catch (const std::exception& e) {

    ExecEnv::log().error("NegativeBinomial::pdf; boost error: {}", e.what());
    return 0.0;

  }

}

double kel::NegativeBinomial::cdf(std::size_t count) const {

  try {

    return bm::cdf(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_), static_cast<double>(count));

  } catch (const std::exception& e) {

    ExecEnv::log().error("NegativeBinomial::cdf; boost error: {}", e.what());
    return 0.0;

  }

}

std::size_t kel::NegativeBinomial::quantile(double p) const {

  try {

    return static_cast<std::size_t>(bm::quantile(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_), p));

  } catch (const std::exception& e) {

    ExecEnv::log().error("NegativeBinomial::quantile; boost error: {}", e.what());
    return 0;

  }

}

double kel::NegativeBinomial::mean() const {

  try {

    return bm::mean(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_));

  } catch (const std::exception& e) {

    ExecEnv::log().error("NegativeBinomial::mean; boost error: {}", e.what());
    return 0.0;

  }

}