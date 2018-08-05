
#ifndef KGD_DISTRIBUTION_H
#define KGD_DISTRIBUTION_H


#include <memory>
#include <random>


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Various random number generators, with the 64 bit Mersenne Twister as the entropy source.
// Note that copy semantics have been disabled as there seems no reasonable reason to do this, and is probably dangerous.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using EntropyGenerator = std::mt19937_64;
class RandomEntropySource
{

public:

  RandomEntropySource() : generator_(rd_()) {}
  RandomEntropySource(const RandomEntropySource&) = delete;
  virtual ~RandomEntropySource() = default;

  RandomEntropySource& operator=(const RandomEntropySource&) = delete;

  EntropyGenerator& generator() const { return generator_; }

private:

  // Order is important! The random device must be created first.
  std::random_device rd_;
  mutable EntropyGenerator generator_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The oxymoronic DeterministicEntropySource object starts the Mersenne Twister from a known point and
// commits the code to a deterministic path. This may be useful for debugging. The start point (seed) can be varied.
// Production code should always use the RandomEntropySource object.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DeterministicEntropySource
{

public:

  DeterministicEntropySource(size_t seed = 1111) : generator_(seed) {}
  DeterministicEntropySource(const DeterministicEntropySource&) = delete;
  virtual ~DeterministicEntropySource() = default;

  DeterministicEntropySource& operator=(const DeterministicEntropySource&) = delete;

  EntropyGenerator& generator() const { return generator_; }

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
// Double floating random on the interval [lower_bound, upper_bound]
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UniformRealDistribution
{

public:

  UniformRealDistribution(double lower_bound, double upper_bound) : uniform_real_(lower_bound, upper_bound) {}
  UniformRealDistribution(const UniformRealDistribution&) = delete;
  virtual ~UniformRealDistribution() = default;

  UniformRealDistribution& operator=(const UniformRealDistribution&) = delete;

  double random(EntropyGenerator &source) const { return uniform_real_(source); }

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
// Integer random numbers between, and including, lower_bound to upper_bound.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class UniformIntegerDistribution
{

public:

  UniformIntegerDistribution(size_t lower_bound, size_t upper_bound) : uniform_integer_(lower_bound, upper_bound) {}
  virtual ~UniformIntegerDistribution() = default;

  size_t random(EntropyGenerator &source) const { return uniform_integer_(source); }

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

  double random(EntropyGenerator &source) const { return normal_real_(source); }

  static double pdf(double x, double mean, double std_dev);

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
// Beta distributed random number implemented as a ratio of gamma random variates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BetaDistribution
{

public:

  BetaDistribution(double a, double b) : x_(a, 1), y_(b, 1) {}
  BetaDistribution(const BetaDistribution&) = delete;
  virtual ~BetaDistribution() = default;

  BetaDistribution& operator=(BetaDistribution&) = delete;

  double random(EntropyGenerator &source) const {  double x = x_(source);  return x / (x + y_(source)); }

  static double logInverseBetaFunction(double a, double b);

  static double logPartialPdf(double x, double a, double b);

  static double logPdf(double x, double a, double b);

  static double pdf(double x, double a, double b);


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

  static double pdf(size_t n, size_t k, double alpha, double beta);
  // Uses loggamma functions.
  static double logPdf(double n, double k, double alpha, double beta);

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

  size_t random(EntropyGenerator &source) const { return binomial_integer_(source); }

  static double pdf(size_t n, size_t k, double prob_success);

private:

  mutable std::binomial_distribution<> binomial_integer_;

};


}   // organization level namespace
}   // project level namespace



#endif
