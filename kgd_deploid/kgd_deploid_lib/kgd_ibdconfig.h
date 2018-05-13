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
namespace deploid {          // project level namespace




// The IBDconfiguration is used for index, which should be non-negative, use int, any thing below zero should throw.
class IBDconfiguration {
#ifdef UNITTEST
  friend class TestIBDconfig;
  friend class TestHprior;
#endif

  friend class Hprior;

  friend class McmcMachinery;

  IBDconfiguration() = default;
  ~IBDconfiguration() = default;

  size_t kStrain_;
  std::vector<std::vector<int> > op;
  std::vector<std::vector<int> > pairToEmission;
  std::vector<std::vector<size_t> > pairList;
  std::vector<std::vector<int> > states;
  std::vector<size_t> effectiveK;

  void buildIBDconfiguration(size_t k = 5);

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  size_t stateSize() const { return states.size(); }

  void enumerateOp();

  void makePairList();

  void makePairToEmission();

  void findUniqueState();

  void findEffectiveK();

  std::vector<int> makeTmpRow();

  std::vector<size_t> findWhichIsOne(std::vector<int> tmpOp);

  std::vector<std::string> getIBDconfigureHeader();

  static int nchoose2(int n);

  static bool twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2);

  static std::vector<std::vector<int> > unique(std::vector<std::vector<int> > &mat);

  static std::vector<int> convertIntToBinary(int x, size_t len);

  static std::vector<std::vector<int> > enumerateBinaryMatrixOfK(size_t k);

};



}   // organization level namespace
}   // project level namespace



#endif
