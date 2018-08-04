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

#ifndef KGD_MERSENNE_TWISTER_H
#define KGD_MERSENNE_TWISTER_H

#include <random>
#include "kgd_random_generator.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class MersenneTwister : public RandomGenerator
{

 public:

  MersenneTwister() : mt_(rd_()) {}
  ~MersenneTwister() = default;

  double unitUniformRand() { return unif_(mt_); }

 private:

  // Order is important! The random device must be created first.
  std::random_device rd_;
  std::mt19937_64 mt_;
  std::uniform_real_distribution<> unif_;

};





}   // organization level namespace
}   // project level namespace



#endif
