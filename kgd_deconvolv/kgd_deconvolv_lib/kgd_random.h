
#ifndef KGD_RANDOM_H
#define KGD_RANDOM_H


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
  ~RandomEntropySource() = default;

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
  ~DeterministicEntropySource() = default;

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


class RandomUniform
{

public:

  RandomUniform(double lower_bound, double upper_bound) : uniform_real_(lower_bound, upper_bound) {}
  RandomUniform(const RandomUniform&) = delete;
  virtual ~RandomUniform() = default;

  RandomUniform& operator=(const RandomUniform&) = delete;

  double generate(EntropyGenerator& source) const { return uniform_real_(source); }

private:

  mutable std::uniform_real_distribution<> uniform_real_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Random numbers on the unit interval [0, 1]
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomUnitUniform : public RandomUniform {

public:

  RandomUnitUniform() : RandomUniform(0.0, 1.0) {}
  ~RandomUnitUniform() override = default;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Integer random numbers between, and including, lower_bound to upper_bound.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomInteger
{

public:

  RandomInteger(size_t lower_bound, size_t upper_bound) : uniform_integer_(lower_bound, upper_bound) {}
  virtual ~RandomInteger() = default;

  size_t generate(EntropyGenerator& source) const { return uniform_integer_(source); }

private:

  mutable std::uniform_int_distribution<size_t> uniform_integer_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Normally distributed random numbers with mean and std_deviation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomNormal
{

public:

  RandomNormal(double mean, double std_deviation) : normal_real_(mean, std_deviation) {}
  RandomNormal(const RandomNormal&) = delete;
  virtual ~RandomNormal() = default;

  RandomNormal& operator=(RandomNormal&) = delete;

  double generate(EntropyGenerator& source) const { return normal_real_(source); }

private:

  mutable std::normal_distribution<> normal_real_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Standard normally distributed random numbers with mean = 0.0 and std_deviation = 1.0.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomStdNormal : public RandomNormal {

public:

  RandomStdNormal() : RandomNormal(0.0, 1.0) {}
  ~RandomStdNormal() override = default;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Beta distributed random number implemented as a ratio of gamma random variates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomBeta
{

public:

  RandomBeta(double a, double b) : x_(a, 1), y_(b, 1) {}
  RandomBeta(const RandomNormal&) = delete;
  virtual ~RandomBeta() = default;

  RandomNormal& operator=(RandomNormal&) = delete;

  double generate(EntropyGenerator& source) const {  double x = x_(source);  return x / (x + y_(source)); }

private:

  mutable std::gamma_distribution<> x_;
  mutable std::gamma_distribution<> y_;

};



}   // organization level namespace
}   // project level namespace



#endif
