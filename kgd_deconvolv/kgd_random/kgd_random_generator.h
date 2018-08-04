/*
 * kgd_deconvolv is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deconvolv.
 *
 * kgd_deconvolv is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef KGD_RANDOM_GENERATOR_H
#define KGD_RANDOM_GENERATOR_H


#include <cassert>
#include <cmath>
#include <memory>
#include <random>


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class RandomGenerator
{

 public:

  RandomGenerator() = default;
  virtual ~RandomGenerator() = default;

  virtual double unitUniformRand() = 0;
  int sampleInt(const int max_value) {  return(static_cast<int>(this->unitUniformRand()*max_value));  }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Various random number generators, with the 64 bit Mersenne Twister as the entropy source.
// Note that copy semantics have been disabled as there seems no reasonable reason to do this, and may be dangerous.
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

  DeterministicEntropySource(size_t seed = 1234) : generator_(seed) {}
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
// Double float normally distributed random numbers with mean and std_deviation.
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
// Double float standard normally distributed random numbers with mean = 0.0 and std_deviation = 1.0.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomStdNormal : public RandomNormal {

public:

  RandomStdNormal() : RandomNormal(0.0, 1.0) {}
  ~RandomStdNormal() override = default;

};



}   // organization level namespace
}   // project level namespace



#endif
