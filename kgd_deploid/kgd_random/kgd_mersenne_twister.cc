/*
 * kgd_deploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deploid.
 *
 * kgd_deploid is free software: you can redistribute it and/or modify
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

#include "kgd_mersenne_twister.h"


namespace kgd = kellerberrin::deploid;


void kgd::MersenneTwister::construct_common(const size_t seed){

  unif_ = std::uniform_real_distribution<>(0, 1);
  set_seed(seed);

}

kgd::MersenneTwister::MersenneTwister() {

  construct_common(generateRandomSeed());

}

kgd::MersenneTwister::MersenneTwister(const size_t seed){

  construct_common(seed);

}

kgd::MersenneTwister::MersenneTwister(const bool use_seed, size_t seed){

  if (!use_seed) seed = generateRandomSeed();
  construct_common(seed);

}

/**
 * @brief Generates a random seed using entropy provided by the operating
 * system.
 *
 * @return A random int between 0 and 2^32
 */
size_t kgd::MersenneTwister::generateRandomSeed() const {

  std::random_device rd;
  std::uniform_int_distribution<size_t> dist(0, 4294967295); // 0 - 2^32-1

  return(dist(rd));

}

void kgd::MersenneTwister::set_seed(const size_t seed) {

  RandomGenerator::set_seed(seed);
  mt_ = std::mt19937_64(seed);
  initializeUnitExponential();

}

