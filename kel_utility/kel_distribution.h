
#ifndef KGD_DISTRIBUTION_H
#define KGD_DISTRIBUTION_H


#include <memory>
#include <random>


namespace kellerberrin {    // organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Random number generators and Pdfs for the Uniform, Discrete Uniform, Guassian, Binomial and BetaBinomial distributions.
// Note that object copy semantics have been disabled as there seems no reasonable reason to do this,
// and in the case of the random number generators, is probably dangerous.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Use the 64 bit Mersenne Twister as the random number entropy source.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using EntropyGenerator = std::mt19937_64;
class RandomEntropySource
{

public:

  RandomEntropySource() : generator_(rd_()) {}
  RandomEntropySource(const RandomEntropySource&) = delete;
  virtual ~RandomEntropySource() = default;

  RandomEntropySource& operator=(const RandomEntropySource&) = delete;

  [[nodiscard]] EntropyGenerator& generator() const { return generator_; }

private:

  // Order is important! The random device must be created first.
  std::random_device rd_;
  mutable EntropyGenerator generator_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The oxymoronic DeterministicEntropySource object starts the Mersenne Twister from a known point and
// commits the code to a deterministic path. This may be useful for debugging. The start point (seed) can be varied.
// Production code should always use the RandomEntropySource object defined above.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DeterministicEntropySource
{

public:

  explicit DeterministicEntropySource(size_t seed = 1111) : generator_(seed) {}
  DeterministicEntropySource(const DeterministicEntropySource&) = delete;
  virtual ~DeterministicEntropySource() = default;

  [[nodiscard]] DeterministicEntropySource& operator=(const DeterministicEntropySource&) = delete;

  [[nodiscard]] EntropyGenerator& generator() const { return generator_; }

private:

  mutable EntropyGenerator generator_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Select the entropy source and recompile.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define RANDOM_ENTROPY_    // Must be enabled for production.
#ifdef RANDOM_ENTROPY_

using EntropySource = RandomEntropySource;

#else

using EntropySource = DeterministicEntropySource;

#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Real (double) random number generator on the interval [lower_bound, upper_bound]
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UniformRealDistribution
{

public:

  UniformRealDistribution(double lower_bound, double upper_bound) : uniform_real_(lower_bound, upper_bound) {}
  UniformRealDistribution(const UniformRealDistribution&) = delete;
  virtual ~UniformRealDistribution() = default;

  UniformRealDistribution& operator=(const UniformRealDistribution&) = delete;

  [[nodiscard]] double random(EntropyGenerator &source) const { return uniform_real_(source); }

private:

  mutable std::uniform_real_distribution<> uniform_real_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Random numbers on the unit interval [0, 1]
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UniformUnitDistribution : public UniformRealDistribution {

public:

  UniformUnitDistribution() : UniformRealDistribution(0.0, 1.0) {}
  virtual ~UniformUnitDistribution() override = default;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Positive integer random numbers between, and including, lower_bound to upper_bound.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UniformIntegerDistribution
{

public:

  UniformIntegerDistribution(size_t lower_bound, size_t upper_bound) : uniform_integer_(lower_bound, upper_bound) {}
  virtual ~UniformIntegerDistribution() = default;

  [[nodiscard]] size_t random(EntropyGenerator &source) const { return uniform_integer_(source); }

private:

  mutable std::uniform_int_distribution<size_t> uniform_integer_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A boolean coin-flip object.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomBoolean
{

public:

  RandomBoolean() : binary_integer_(0, 1) {}
  virtual ~RandomBoolean() = default;

  [[nodiscard]] bool random(EntropyGenerator &source) const { return binary_integer_.random(source) > 0; }

private:

  UniformIntegerDistribution binary_integer_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Normally distributed random numbers with mean and std_deviation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class NormalDistribution
{

public:

  NormalDistribution(double mean, double std_deviation) : normal_real_(mean, std_deviation) {}
  NormalDistribution(const NormalDistribution&) = delete;
  virtual ~NormalDistribution() = default;

  NormalDistribution& operator=(NormalDistribution&) = delete;

  [[nodiscard]] double random(EntropyGenerator &source) const { return normal_real_(source); }

  [[nodiscard]] static double pdf(double x, double mean, double std_dev);

private:

  mutable std::normal_distribution<> normal_real_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard normally distributed random numbers with mean = 0.0 and std_deviation = 1.0.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class StdNormalDistribution : public NormalDistribution {

public:

  StdNormalDistribution() : NormalDistribution(0.0, 1.0) {}
  virtual ~StdNormalDistribution() override = default;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gamma distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GammaDistribution
{

public:

  GammaDistribution(double a, double b) : x_(a, b), _a(a), _b(b) {}
  GammaDistribution(const GammaDistribution&) = delete;
  virtual ~GammaDistribution() = default;

  GammaDistribution& operator=(GammaDistribution&) = delete;

  [[nodiscard]] double random(EntropyGenerator &source) const {  return x_(source); }

  [[nodiscard]] double pdf(double x) const;

  [[nodiscard]] double cdf(double x) const;

  [[nodiscard]] double quantile(double p) const;

private:

  mutable std::gamma_distribution<> x_;
  double _a;
  double _b;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta distributed random numbers implemented as a ratio of gamma random variates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BetaDistribution
{

public:

  BetaDistribution(double a, double b) : x_(a, 1), y_(b, 1) {}
  BetaDistribution(const BetaDistribution&) = delete;
  virtual ~BetaDistribution() = default;

  BetaDistribution& operator=(BetaDistribution&) = delete;

  [[nodiscard]] double random(EntropyGenerator &source) const {  double x = x_(source);  return x / (x + y_(source)); }

  [[nodiscard]] static double logInverseBetaFunction(double a, double b);

  [[nodiscard]] static double logPartialPdf(double x, double a, double b);

  [[nodiscard]] static double logPdf(double x, double a, double b);

  [[nodiscard]] static double pdf(double x, double a, double b);

  [[nodiscard]] static double mean(double a, double b);

  [[nodiscard]] static double var(double a, double b);

  [[nodiscard]] static double mode(double a, double b);

private:

  mutable std::gamma_distribution<> x_;
  mutable std::gamma_distribution<> y_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta Binomial Distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BetaBinomialDistribution {

public:

  BetaBinomialDistribution() = delete;
  ~BetaBinomialDistribution() = delete;

  // Uses boost beta functions.
  [[nodiscard]] static double pdf(size_t n, size_t k, double alpha, double beta);
  // Uses boost lgamma functions.
  [[nodiscard]] static double logPdf(double n, double k, double alpha, double beta);
  // No binomial coefficient term.
  [[nodiscard]] static double logPartialPdf(double n, double k, double alpha, double beta);
  // No binomial coefficient term.
  [[nodiscard]] static double partialPdf(double n, double k, double alpha, double beta);

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Binomial Distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BinomialDistribution {

public:

  BinomialDistribution(size_t trials, double prob_success) :  binomial_integer_(trials, prob_success) {}
  BinomialDistribution(const BinomialDistribution &) = delete;
  ~BinomialDistribution() = default;

  BinomialDistribution &operator=(BinomialDistribution &) = delete;

  [[nodiscard]] size_t random(EntropyGenerator &source) const { return binomial_integer_(source); }

  [[nodiscard]] static double pdf(size_t n, size_t k, double prob_success);

  [[nodiscard]] static double cdf(size_t n, double k, double prob_success);

  [[nodiscard]] static double mean(size_t n, double prob_success) { return static_cast<double>(n) * prob_success; }

private:

  mutable std::binomial_distribution<> binomial_integer_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Bernoulli Distribution
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BernoulliDistribution {

public:

  BernoulliDistribution(double prob_success) :  bernoulli_integer_(prob_success) {}
  BernoulliDistribution(const BinomialDistribution &) = delete;
  ~BernoulliDistribution() = default;

  BinomialDistribution &operator=(BinomialDistribution &) = delete;

  [[nodiscard]] bool random(EntropyGenerator &source) const { return bernoulli_integer_(source); }

  [[nodiscard]] static double mean(double prob_success) { return prob_success; }

private:

  mutable std::bernoulli_distribution bernoulli_integer_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hypergeometric distribution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HypergeometricDistribution
{

public:

  HypergeometricDistribution(size_t pop_successes_K, size_t sample_size_n, size_t population_N);
  HypergeometricDistribution(const HypergeometricDistribution&) = delete;
  ~HypergeometricDistribution() = default;

  HypergeometricDistribution& operator=(HypergeometricDistribution&) = delete;

  [[nodiscard]] double pdf(size_t successes_k) const;

  [[nodiscard]] double cdf(size_t successes_k) const;

  [[nodiscard]] double quantile(size_t successes_k) const;

//  The hypergeometric tests below use the hypergeometric distribution to measure the statistical significance
//  of having drawn a sample consisting of a specific number of k successes n total draws (without replacement)
//  from a population of size N containing K successes.

//  The test for over-representation of successes in the sample, the hypergeometric p-value is calculated
//  as the probability of randomly drawing k or more successes from the population in n total draws.
  [[nodiscard]] double upperSingleTailTest(size_t successes_k) const;

//  The test for under-representation, the p-value is the probability of randomly drawing k or fewer successes.
  [[nodiscard]] double lowerSingleTailTest(size_t successes_k) const;

  // Bounds for the number of successes_k in a drawn sample_size_n (without replacement)
  // The number of successes drawn only has support; k in { lowerSuccesses_k, ... , upperSuccesses_k}
  [[nodiscard]] size_t upperSuccesses_k() const { return std::min(pop_successes_K_, sample_size_n_); }
  [[nodiscard]] size_t lowerSuccesses_k() const {

    int64_t lower_success = static_cast<int64_t>(sample_size_n_ + pop_successes_K_) - static_cast<int64_t>(population_N_);
    return static_cast<size_t>(std::max<int64_t>(0, lower_success));

  }

private:

  size_t pop_successes_K_;
  size_t sample_size_n_;
  size_t population_N_;


};


}   // end namespace


#endif
