//
// Created by kellerberrin on 5/08/18.
//

#include "kel_exec_env.h"
#include "kel_distribution.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>


namespace bm = boost::math;
namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Distribution functions use boost special functions.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////


double kel::NormalDistribution::pdf(double x, double mean, double std_dev) {

  static const double inv_sqrt_2pi = 0.3989422804014327;

  double a = (x - mean) / std_dev;

  return (inv_sqrt_2pi / std_dev) * std::exp(-0.5 * a * a);

}


double kel::GammaDistribution::pdf(double x) const {

  return bm::pdf(bm::gamma_distribution<>(_a, _b), x);

}

double kel::GammaDistribution::cdf(double x) const {

  return bm::cdf(bm::gamma_distribution<>(_a, _b), x);

}

double kel::GammaDistribution::quantile(double p) const {

  return bm::quantile(bm::gamma_distribution<>(_a, _b), p);

}



double kel::BetaDistribution::logInverseBetaFunction(double a, double b) {

  assert(a > 0.0);
  assert(b > 0.0);

  double log_gamma = bm::lgamma<double>(a + b) - bm::lgamma<double>(a) - bm::lgamma<double>(b);

  return log_gamma;

}


double kel::BetaDistribution::logPartialPdf(double x, double a, double b) {

  assert(x >= 0 && x <= 1);
  assert(a >= 0);
  assert(b >= 0);

  double ret =  (b - 1) * std::log(1 - x) + (a - 1) * std::log(x);

  return ret;

}


double kel::BetaDistribution::logPdf(double x, double a, double b) {

  assert(x >= 0 && x <= 1);
  assert(a > 0);
  assert(b > 0);

  double ret = logInverseBetaFunction(a, b) + logPartialPdf(x, a, b);

  return ret;

}

double kel::BetaDistribution::pdf(double x, double a, double b) {

  assert(x >= 0 && x <= 1);
  assert(a > 0);
  assert(b > 0);

  double p = bm::tgamma<double>(a + b) / (bm::tgamma<double>(a) * bm::tgamma<double>(b));
  double q = std::pow(1 - x, b - 1) * std::pow(x, a - 1);
  return p * q;

}


double kel::BetaDistribution::mean(double a, double b) {

  assert(a > 0);
  assert(b > 0);

  return a / (a + b);

}


double kel::BetaDistribution::var(double a, double b) {

  assert(a > 0);
  assert(b > 0);

  return (a * b) / ((a + b) * (a + b) * (a + b + 1.0));

}


double kel::BetaDistribution::mode(double a, double b) {

  assert(a > 1);
  assert(b > 1);

  return (a - 1.0) / (a + b - 2.0);

}

double kel::BetaBinomialDistribution::pdf(size_t n, size_t k, double alpha, double beta) {

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


double kel::BetaBinomialDistribution::partialPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0);
  assert(beta > 0);

  double a1 = k + alpha;
  double b1 = n - k + beta;

  double p = bm::beta<double>(a1, b1);
  double q = bm::beta<double>(alpha, beta);

  return (p / q);

}


double kel::BetaBinomialDistribution::logPartialPdf(double n, double k, double alpha, double beta) {

  assert(k <= n);
  assert(alpha > 0.0);
  assert(beta > 0.0);

  double r = n - k;
  double a1 = k + alpha;
  double b1 = r + beta;

  double beta_n = bm::lgamma<double>(a1) + bm::lgamma<double>(b1) - bm::lgamma<double>(a1 + b1);
  double beta_d = bm::lgamma<double>(alpha) + bm::lgamma<double>(beta) - bm::lgamma<double>(alpha + beta);
  double prob = beta_n - beta_d;

  return prob;

}


double kel::BetaBinomialDistribution::logPdf(double n, double k, double alpha, double beta) {

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


// .first is alpha, .second is beta. The raw moments are calculated from the observations and used to calculate alpha and beta.
[[nodiscard]] std::pair<double, double> kel::BetaBinomialDistribution::methodOfMoments(const std::vector<size_t>& observations, size_t n_trials) {

  // calculate the first moment
  size_t obs_sum = std::accumulate(observations.begin(), observations.end(), static_cast<size_t>(0));
  const double m1 = static_cast<double>(obs_sum) / static_cast<double>(observations.size());

  // calculate the 2nd moment
  auto sqr_lambda = [](size_t accumulate, size_t obs)->size_t { return accumulate + (obs * obs); };
  size_t sqr_sum = std::accumulate(observations.begin(), observations.end(), static_cast<size_t>(0), sqr_lambda);
  const double m2 = static_cast<double>(sqr_sum) / static_cast<double>(observations.size());

  const double n = static_cast<double>(n_trials);

  const double a_numer = (n * m1) - m2;
  const double ab_denom = n * ((m2 / m1) - m1 - 1) + m1;
  const double a = a_numer / ab_denom;

  const double b_numer = (n - m1) * (n - (m2/m1));
  const double b = b_numer / ab_denom;

  return {a, b};

}



double kel::BinomialDistribution::pdf(size_t n, size_t k, double prob_success) {

  assert(k <= n and k >= 0);
  assert(prob_success >= 0 and prob_success <= 1.0);

  double coeff = bm::binomial_coefficient<double>(n, k);

  double p = std::pow(prob_success, static_cast<double>(k));

  double q = std::pow((1.0 - prob_success), static_cast<double>(n - k));

  return coeff * p * q;

}


double kel::BinomialDistribution::cdf(size_t n, double k, double prob_success) {

  assert(k <= n and k >= 0);
  assert(prob_success >= 0 and prob_success <= 1.0);

  size_t integer_k = std::floor(k);
  double sum_pdf = 0.0;

  for (size_t index = 0; index <= integer_k; ++index) {

    sum_pdf += cdf(n, index, prob_success);

  }

  return sum_pdf;

}


kel::HypergeometricDistribution::HypergeometricDistribution(size_t pop_successes_K, size_t sample_size_n, size_t population_N) {


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

double kel::HypergeometricDistribution::pdf(size_t successes_k) const {

  if (successes_k > upperSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution::pdf; r_successes k:{} exceeds upper limit :{}",
                         successes_k, upperSuccesses_k());
    successes_k = upperSuccesses_k();

  }

  if (successes_k < lowerSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution::pdf; r_successes k:{} below lower limit :{}",
                         successes_k, lowerSuccesses_k());
    successes_k = lowerSuccesses_k();

  }

  bm::hypergeometric_distribution hypergeometric(pop_successes_K_, sample_size_n_, population_N_);

  return bm::pdf(hypergeometric, successes_k);

}

double kel::HypergeometricDistribution::cdf(size_t successes_k) const {

  if (successes_k > upperSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution(N:{},K:{},n:{},k:{})::cdf; r_successes k exceeds upper limit :{}",
                        population_N_, pop_successes_K_, sample_size_n_, successes_k, upperSuccesses_k());
    successes_k = upperSuccesses_k();

  }

  if (successes_k < lowerSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution(N:{},K:{},n:{},k:{})::cdf; r_successes k exceeds upper limit :{}",
                        population_N_, pop_successes_K_, sample_size_n_, successes_k, lowerSuccesses_k());
    successes_k = lowerSuccesses_k();

  }

  bm::hypergeometric_distribution hypergeometric(pop_successes_K_, sample_size_n_, population_N_);

  return bm::cdf(hypergeometric, successes_k);

}

double kel::HypergeometricDistribution::quantile(size_t successes_k) const {

  if (successes_k > upperSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution::quantile; r_successes k:{} exceeds upper limit :{}",
                         successes_k, upperSuccesses_k());
    successes_k = upperSuccesses_k();

  }

  if (successes_k < lowerSuccesses_k()) {

    ExecEnv::log().warn("HypergeometricDistribution::quantile; r_successes k:{} below lower limit :{}",
                         successes_k, lowerSuccesses_k());
    successes_k = lowerSuccesses_k();

  }

  bm::hypergeometric_distribution hypergeometric(pop_successes_K_, sample_size_n_, population_N_);

  return bm::quantile(hypergeometric, successes_k);

}


size_t kel::HypergeometricDistribution::lowerSuccesses_k() const {

  int64_t lower_success = static_cast<int64_t>(sample_size_n_ + pop_successes_K_) - static_cast<int64_t>(population_N_);
  return static_cast<size_t>(std::max<int64_t>(0, lower_success));

}


double kel::HypergeometricDistribution::upperSingleTailTest(size_t test_value_k) const {

  if (test_value_k > 0) {

    test_value_k = test_value_k - 1;

  }

  return 1.0 - cdf(test_value_k);

}

double kel::HypergeometricDistribution::lowerSingleTailTest(size_t test_value_k) const {

  return cdf(test_value_k);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Poisson distribution. Uses boost for implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::Poisson::pdf(size_t count) const {

  return bm::pdf(bm::poisson_distribution<>(lambda_), count);

}

double kel::Poisson::cdf(size_t count) const  {

  return bm::cdf(bm::poisson_distribution<>(lambda_), count);

}

size_t kel::Poisson::quantile(double quantile) const {

  return bm::quantile(bm::poisson_distribution<>(lambda_), quantile);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Negative Binomial distribution. Uses boost for implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double kel::NegativeBinomial::pdf(size_t count) const {

  return bm::pdf(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_), count);

}

double kel::NegativeBinomial::cdf(size_t count) const {

  return bm::cdf(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_), count);

}

size_t kel::NegativeBinomial::quantile(double quantile) const {

  return bm::quantile(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_), quantile);

}

double kel::NegativeBinomial::mean() const {

  return bm::mean(bm::negative_binomial_distribution<>(r_successes_, p_prob_success_));

}
