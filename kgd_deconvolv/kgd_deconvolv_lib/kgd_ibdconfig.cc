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

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "kgd_deconvolv_app.h"
#include "kgd_combinatorial.h"
#include "kgd_ibdconfig.h"


namespace kgd = kellerberrin::deconvolv;

//#define IBDCONFIG_DEBUG 1

void kgd::IBDconfiguration::buildIBDconfiguration(size_t k) {

  setKstrain(k);

  std::vector< std::vector<size_t>> pairs_list;
  makePairList(pairs_list);

  std::vector< std::vector<size_t>> pairs_to_emission;
  makePairToEmission(pairs_list, pairs_to_emission);

  /// Only unique pairs to emission.
  states_ = Utility::uniqueMatrixColumns(pairs_to_emission);

  findEffectiveK();

#ifdef IBDCONFIG_DEBUG

  ExecEnv::log().info("pairs_to_emission.size:{}, states_ size: {}, states_[0].size:{}, effective_k.size: {}",
                      pairs_to_emission.size(), states_.size(), states_[0].size(), effectiveK_.size());

  size_t index = 0;
  for (auto state: getIBDconfigureHeader()) {

    ExecEnv::log().info("State[{}]: {}, EffectiveK: {}", index, state, effectiveK_[index]);
    index++;

  }

#endif

}


void kgd::IBDconfiguration::makePairList(std::vector< std::vector<size_t>>& pairs_list) {
  //#Make a map of pairs to pair value
  /// Creates a pair list (a,b) where we have k items. This will be (k * (k -1 ) /2) paired items.

  assert(pairs_list.empty());

  for (size_t i = 0; i < kStrain(); i++) { // 0-indexed

    for (size_t j = i + 1; j < kStrain(); j++) {

      pairs_list.push_back(std::vector<size_t>({i, j}));

    }

  }

  assert(pairs_list.size() == (size_t) nchoose2(kStrain()));

}


void kgd::IBDconfiguration::makePairToEmission(const std::vector< std::vector<size_t>>& pairs_list,
                                               std::vector< std::vector<size_t>>& pairs_to_emission) {

  assert(pairs_to_emission.empty());

  /// Generates a matrix of 2^ (k * (k -1 ) /2) rows where each row is a binary representation of all [0, 2^ (k * (k -1 ) /2))
  BinaryPermutations permutations(nchoose2(kStrain()));

  /// For all permuted pair combinations.
  for (size_t bit_pattern = 0; bit_pattern < permutations.bitArraySize(); ++bit_pattern) {

    /// Create a vector of size k which is pre-initialized to [0, ..., (k -1)]
    std::vector<size_t> enumerated_array = makeEnumeratedArray();

    /// For each permutation create an array of indexes of active pairs.
    std::vector<size_t> reverse_bit_array = permutations.reverseBitArray(bit_pattern);

    std::vector<size_t> active_pairs_array = activePairsArray(reverse_bit_array);

    /// If we have active pairs, create an active pair list.
    if (active_pairs_array.size() > 0) {

      std::vector<std::vector<size_t> > active_pair_list;

      for (size_t j = 0; j < active_pairs_array.size(); j++) {

        active_pair_list.push_back(pairs_list[active_pairs_array[j]]);

      }

      int tmpIndex = (active_pair_list.size() - 1);
      /// For all active pairs (a, b), array[a] = array[b]
      while (tmpIndex >= 0) {

        enumerated_array[active_pair_list[tmpIndex][0]] = enumerated_array[active_pair_list[tmpIndex][1]];
        tmpIndex--;

      }

    }

    pairs_to_emission.push_back(enumerated_array);

  }

}


std::vector<size_t> kgd::IBDconfiguration::makeEnumeratedArray() {

  std::vector<size_t> ret(kStrain());

  for (size_t i = 0; i < ret.size(); i++) {

    ret[i] = i;

  }

  return ret;

}


std::vector<size_t> kgd::IBDconfiguration::activePairsArray(const std::vector<size_t>& pair_permute_row) {

/// This function returns an array of indexes of active pairs for a pair_permute_row.
  std::vector<size_t> ret;

  for (size_t i = 0; i < pair_permute_row.size(); i++) {

    if (pair_permute_row[i] == 1) {

      ret.push_back(i);

    }

  }

  return ret;

}


void kgd::IBDconfiguration::findEffectiveK() {

  assert(effectiveK_.empty());

  for (auto state : states_) {

    std::set<size_t> tmpSet(state.begin(), state.end());

    effectiveK_.push_back(tmpSet.size());

  }

  assert(effectiveK_.size() == states_.size());

}


std::vector<std::string> kgd::IBDconfiguration::getIBDconfigureHeader() const {

  std::vector<std::string> ret;

  for (size_t i = 0; i < states_.size(); i++) {

    std::string tmp;

    for (size_t j = 0; j < states_[i].size(); j++) {

      std::stringstream tmp_ss;
      tmp_ss << states_[i][j];
      tmp += tmp_ss.str() + ((j < (states_[i].size() - 1)) ? "-" : "");

    }

    ret.push_back(tmp);

  }

  return ret;

}


size_t kgd::IBDconfiguration::nchoose2(size_t n) {

  if (n < 2) {

    DeconvolvApp::log().critical("IBDconfiguration::nchoose2(n); Domain error n: {} must be >= 2", n);

  }

  size_t ret = n * (n - 1) / 2;

  return ret;

}

