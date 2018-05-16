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



#ifndef KGD_IBD_H
#define KGD_IBD_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_mersenne_twister.h"
#include "kgd_deploid_io.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


// The IBDconfiguration is used for index, which should be non-negative, use int, any thing below zero should throw.
class IBDconfiguration {

#ifdef UNITTEST
  friend class TestIBDconfig;
  friend class TestHprior;
#endif

public:

  IBDconfiguration() = default;
  ~IBDconfiguration() = default;

  void buildIBDconfiguration(size_t k = 5);

  static std::vector<std::vector<int> > unique(std::vector<std::vector<int> > &mat);
  static std::vector<std::vector<int> > enumerateBinaryMatrixOfK(size_t k);

  const std::vector<size_t>& effectiveK() const { return effectiveK_; }
  const std::vector<std::vector<int> >& states() const { return states_; }
  std::vector<std::string> getIBDconfigureHeader() const;

private:

  size_t kStrain_;
  std::vector<std::vector<int> > pairs_permutation_;
  std::vector<std::vector<int> > pair_to_emission_;
  std::vector<std::vector<size_t> > pair_list_;
  std::vector<std::vector<int> > states_;
  std::vector<size_t> effectiveK_;


  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  size_t stateSize() const { return states_.size(); }

  void enumerateOp();

  void makePairList();

  void makePairToEmission();

  void findUniqueState();

  void findEffectiveK();

  std::vector<int> makeEnumeratedArray();

  std::vector<size_t> activePairsArray(std::vector<int> pair_permute_row);

  static int nchoose2(int n);

  static bool twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2);

  static std::vector<int> convertIntToBinary(int x, size_t len);


};



}   // organization level namespace
}   // project level namespace



#endif
