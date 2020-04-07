
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
// Beta distributed random numbers implemented as a ratio of gamma random variates.
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
  virtual ~BinomialDistribution() = default;

  BinomialDistribution &operator=(BinomialDistribution &) = delete;

  [[nodiscard]] size_t random(EntropyGenerator &source) const { return binomial_integer_(source); }

  [[nodiscard]] static double pdf(size_t n, size_t k, double prob_success);

private:

  mutable std::binomial_distribution<> binomial_integer_;

};


}   // end namespace


#endif
