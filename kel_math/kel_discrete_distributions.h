//
// Created by kellerberrin on 10/8/26.
//

#ifndef KEL_DISCRETE_DISTRIBUTIONS_H
#define KEL_DISCRETE_DISTRIBUTIONS_H


#include <cstddef>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include "kel_entropy_source.h"


namespace kellerberrin {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Discrete distributions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Positive integer random numbers between, and including, lower_bound to upper_bound.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UniformIntegerDistribution final {
public:

  explicit UniformIntegerDistribution(std::size_t lower_bound, std::size_t upper_bound)
        : dist_(lower_bound, upper_bound) {}

  UniformIntegerDistribution(const UniformIntegerDistribution&) = delete;
  UniformIntegerDistribution(UniformIntegerDistribution&&) = default;
  UniformIntegerDistribution& operator=(const UniformIntegerDistribution&) = delete;
  UniformIntegerDistribution& operator=(UniformIntegerDistribution&&) = default;
  ~UniformIntegerDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] std::size_t random(URBG& source) const { return dist_(source); }

private:

  mutable std::uniform_int_distribution<std::size_t> dist_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A boolean coin-flip object.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RandomBoolean final {
public:

  RandomBoolean() = default;

  RandomBoolean(const RandomBoolean&) = delete;
  RandomBoolean(RandomBoolean&&) = default;
  RandomBoolean& operator=(const RandomBoolean&) = delete;
  RandomBoolean& operator=(RandomBoolean&&) = default;
  ~RandomBoolean() = default;

  template<std::uniform_random_bit_generator URBG>
    [[nodiscard]] bool random(URBG& source) const { return dist_(source); }

private:

  mutable std::bernoulli_distribution dist_{0.5};

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Binomial Distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BinomialDistribution final {
public:
  explicit BinomialDistribution(std::size_t trials, double prob_success)
        : dist_(trials, prob_success) {}

  BinomialDistribution(const BinomialDistribution&) = delete;
  BinomialDistribution(BinomialDistribution&&) = default;
  BinomialDistribution& operator=(const BinomialDistribution&) = delete;
  BinomialDistribution& operator=(BinomialDistribution&&) = default;
  ~BinomialDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] std::size_t random(URBG& source) const { return dist_(source); }

  [[nodiscard]] static double pdf(std::size_t n, std::size_t k, double prob_success);
  [[nodiscard]] static double cdf(std::size_t n, double k, double prob_success);
  [[nodiscard]] static double mean(std::size_t n, double prob_success) {
        return static_cast<double>(n) * prob_success;
  }

private:

  mutable std::binomial_distribution<std::size_t> dist_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Bernoulli Distribution
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BernoulliDistribution final {
public:

  explicit BernoulliDistribution(double prob_success) : dist_(prob_success) {}

  BernoulliDistribution(const BernoulliDistribution&) = delete;
  BernoulliDistribution(BernoulliDistribution&&) = default;
  BernoulliDistribution& operator=(const BernoulliDistribution&) = delete;
  BernoulliDistribution& operator=(BernoulliDistribution&&) = default;
  ~BernoulliDistribution() = default;

  template<std::uniform_random_bit_generator URBG>
  [[nodiscard]] bool random(URBG& source) const { return dist_(source); }

  [[nodiscard]] static double mean(double prob_success) { return prob_success; }

private:

  mutable std::bernoulli_distribution dist_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta Binomial Distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Static-only helper: a namespace is cleaner than a class with deleted ctor/dtor.
namespace BetaBinomialDistribution {

// Uses boost beta functions.
  [[nodiscard]] double pdf(std::size_t n, std::size_t k, double alpha, double beta);
// Uses boost lgamma functions.
  [[nodiscard]] double logPdf(double n, double k, double alpha, double beta);
// No binomial coefficient term.
  [[nodiscard]] double logPartialPdf(double n, double k, double alpha, double beta);
// No binomial coefficient term.
  [[nodiscard]] double partialPdf(double n, double k, double alpha, double beta);
// .first is alpha, .second is beta. The raw moments are calculated from the observations and used to calculate alpha and beta.
  [[nodiscard]] std::pair<double, double> methodOfMoments( const std::vector<std::size_t>& observations, std::size_t n_trials);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distribution objects that only provide pdf/cdf/quantile
// (no internal RNG state, so copying is safe and useful).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hypergeometric distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HypergeometricDistribution final {
public:
  explicit HypergeometricDistribution(std::size_t pop_successes_K,
                                      std::size_t sample_size_n,
                                      std::size_t population_N);

  HypergeometricDistribution(const HypergeometricDistribution&) = default;
  HypergeometricDistribution(HypergeometricDistribution&&) = default;
  HypergeometricDistribution& operator=(const HypergeometricDistribution&) = default;
  HypergeometricDistribution& operator=(HypergeometricDistribution&&) = default;
  ~HypergeometricDistribution() = default;

  [[nodiscard]] double pdf(std::size_t successes_k) const;
  [[nodiscard]] double cdf(std::size_t successes_k) const;
  [[nodiscard]] double quantile(std::size_t successes_k) const;

  //  The hypergeometric tests below use the hypergeometric distribution to measure the statistical significance
  //  of having drawn a sample consisting of a specific number of k r_successes n total draws (without replacement)
  //  from a population of size N containing K r_successes.

  [[nodiscard]] double upperSingleTailTest(std::size_t successes_k) const;
  //  The test for under-representation, the p-value is the prob failure of randomly drawing k or fewer r_successes.
  [[nodiscard]] double lowerSingleTailTest(std::size_t successes_k) const;

  // Bounds for the number of successes_k in a drawn sample_size_n (without replacement)
  // The number of r_successes drawn only has support; k in { lowerSuccesses_k, ... , upperSuccesses_k}
  [[nodiscard]] std::size_t upperSuccesses_k() const noexcept { return std::min(pop_successes_K_, sample_size_n_);}
  [[nodiscard]] std::size_t lowerSuccesses_k() const;

private:

  std::size_t pop_successes_K_;
  std::size_t sample_size_n_;
  std::size_t population_N_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Poisson distribution. Uses boost for implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Poisson final {
public:

  explicit Poisson(double lambda) : lambda_(lambda) {}

  Poisson(const Poisson&) = default;
  Poisson(Poisson&&) = default;
  Poisson& operator=(const Poisson&) = default;
  Poisson& operator=(Poisson&&) = default;
  ~Poisson() = default;

  [[nodiscard]] double pdf(std::size_t count) const;
  [[nodiscard]] double cdf(std::size_t count) const;
  [[nodiscard]] std::size_t quantile(double p) const;
  [[nodiscard]] double mean() const noexcept { return lambda_; }

private:

  double lambda_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The Negative Binomial distribution. Uses boost for implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class NegativeBinomial final {
public:
  explicit NegativeBinomial(double r_successes, double p_prob_success)
        : r_successes_(r_successes), p_prob_success_(p_prob_success) {}

  NegativeBinomial(const NegativeBinomial&) = default;
  NegativeBinomial(NegativeBinomial&&) = default;
  NegativeBinomial& operator=(const NegativeBinomial&) = default;
  NegativeBinomial& operator=(NegativeBinomial&&) = default;
  ~NegativeBinomial() = default;

  [[nodiscard]] double pdf(std::size_t count) const;
  [[nodiscard]] double cdf(std::size_t count) const;
  [[nodiscard]] std::size_t quantile(double p) const;
  [[nodiscard]] double mean() const;

  [[nodiscard]] double r_successes() const noexcept { return r_successes_; }
  [[nodiscard]] double p_probsuccess() const noexcept { return p_prob_success_; }

private:

  double r_successes_;
  double p_prob_success_;

};

}   // namespace kellerberrin


#endif //KEL_DISCRETE_DISTRIBUTIONS_H
