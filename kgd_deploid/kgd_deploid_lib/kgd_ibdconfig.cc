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

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <kgl_exec_env.h>
#include "kgd_ibdconfig.h"


namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;



void kgd::IBDconfiguration::buildIBDconfiguration(size_t k) {

  setKstrain(k);
  enumerateOp();
  makePairList();
  makePairToEmission();
  findUniqueState();
  findEffectiveK();

}



void kgd::IBDconfiguration::enumerateOp() {
  //#For each configuration, identify which pairs are IBD

  /// Generates a matrix of 2^ (k * (k -1 ) /2) rows where each row is a binary representation of all [0, 2^ (k * (k -1 ) /2))
  /// This causes an exponential tial explosition.
  pairs_permutation_ = enumerateBinaryMatrixOfK(nchoose2(kStrain()));

}


void kgd::IBDconfiguration::makePairList() {
  //#Make a map of pairs to pair value
  /// Creates a pair list (a,b) where we have k items. This will be (k * (k -1 ) /2) paired items.

  assert(pair_list_.size() == 0);

  for (size_t i = 0; i < kStrain(); i++) { // 0-indexed

    for (size_t j = i + 1; j < kStrain(); j++) {

      pair_list_.push_back(std::vector<size_t>({i, j}));

    }

  }

  assert(pair_list_.size() == (size_t) nchoose2(kStrain()));

}


void kgd::IBDconfiguration::makePairToEmission() {

  assert(pair_to_emission_.size() == 0);

  /// For all permuted pair combinations.
  for (auto pair_permute_row : pairs_permutation_) {

    /// Create a vector of size k which is pre-initialized to [0, ..., (k -1)]
    std::vector<int> enumerated_array = makeEnumeratedArray();

    /// For each permutation create an array of indexes of active pairs.
    std::vector<size_t> active_pairs_array = activePairsArray(pair_permute_row);

    /// If we have active pairs, create an active pair list.
    if (active_pairs_array.size() > 0) {

      std::vector<std::vector<size_t> > active_pair_list;

      for (size_t j = 0; j < active_pairs_array.size(); j++) {

        active_pair_list.push_back(pair_list_[active_pairs_array[j]]);

      }

      int tmpIndex = (active_pair_list.size() - 1);
      /// For all active pairs (a, b), array[a] = array[b]
      while (tmpIndex >= 0) {

        enumerated_array[active_pair_list[tmpIndex][0]] = enumerated_array[active_pair_list[tmpIndex][1]];
        tmpIndex--;

      }

    }

    pair_to_emission_.push_back(enumerated_array);

  }

}


std::vector<int> kgd::IBDconfiguration::makeEnumeratedArray() {

  std::vector<int> ret(kStrain());

  for (size_t i = 0; i < ret.size(); i++) {

    ret[i] = (int) i;

  }

  return ret;

}


std::vector<size_t> kgd::IBDconfiguration::activePairsArray(std::vector<int> pair_permute_row) {

/// This function returns an array of indexes of active pairs for a pair_permute_row.
  std::vector<size_t> ret;

  for (size_t i = 0; i < pair_permute_row.size(); i++) {

    if (pair_permute_row[i] == 1) {

      ret.push_back(i);

    }

  }

  return ret;

}


void kgd::IBDconfiguration::findUniqueState() {

  assert (states_.size() == 0);

  /// Only unique pairs to to emission.
  states_ = unique(pair_to_emission_);

  kgl::ExecEnv::log().info("pair_to_emission_ size:{}, states_ size: {}, states_[0] size:{}", pair_to_emission_.size(), states_.size(), states_[0].size());

}


void kgd::IBDconfiguration::findEffectiveK() {

  assert(effectiveK_.size() == 0);

  for (auto state : states_) {

    std::set<int> tmpSet(state.begin(), state.end());

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


std::vector<std::vector<int> > kgd::IBDconfiguration::unique(std::vector<std::vector<int> > &mat) {

  std::vector<std::vector<int> > ret;

  ret.push_back(mat[0]);

  for (size_t i = 1; i < mat.size(); i++) {

    bool aNewState = true;

    for (auto state : ret) {

      if (twoVectorsAreSame(state, mat[i])) {

        aNewState = false;
        break;

      }

    }

    if (aNewState) {

      ret.push_back(mat[i]);

    }

  }

  return ret;

}


std::vector<std::vector<int> > kgd::IBDconfiguration::enumerateBinaryMatrixOfK(size_t k) {
  // This function enumerate all possible binary combinations of k elements

  //TODO: Exponential explosion causes memory exhaustion for (k > 5). Replace/remove this logic.
  int ksq = pow(2, k);

  std::vector<std::vector<int> > ret;

  for (int i = 0; i < ksq; i++) {

    ret.push_back(convertIntToBinary(i, k));

  }

  return ret;

}


std::vector<int> kgd::IBDconfiguration::convertIntToBinary(int x, size_t len) {

  std::vector<int> ret(len);

  size_t idx = 0;

  while (x) {

    ret[idx] = (x & 1) ? 1 : 0;
    idx++;

    if (idx > len) {

      throw OutOfVectorSize();
    }

    x >>= 1;

  }

  std::reverse(ret.begin(), ret.end());

  return ret;
}


int kgd::IBDconfiguration::nchoose2(int n) {

  if (n < 2) {

    throw InvalidInput("Input must be at least 2!");

  }

  int ret = n * (n - 1) / 2;
  return ret;

}


bool kgd::IBDconfiguration::twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2) {

  if (vec1.size() != vec2.size()) {

    throw InvalidInput("Input vectors have different length!");

  }

  bool ret = true;

  for (size_t i = 0; i < vec1.size(); i++) {

    if (vec1[i] != vec2[i]) {

      ret = false;
      break;

    }

  }

  return ret;

}

