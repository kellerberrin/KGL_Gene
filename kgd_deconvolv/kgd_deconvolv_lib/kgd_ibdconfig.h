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
#include "kgd_deploid_io.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


// The IBDconfiguration is used for indexing.

class IBDconfiguration {

public:

  IBDconfiguration() = default;
  ~IBDconfiguration() = default;

  void buildIBDconfiguration(size_t k = 5);

  const std::vector<size_t>& effectiveK() const { return effectiveK_; }
  const std::vector<std::vector<size_t> >& states() const { return states_; }
  std::vector<std::string> getIBDconfigureHeader() const;

private:

  size_t kStrain_;
  std::vector<std::vector<size_t> > states_;
  std::vector<size_t> effectiveK_;


  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  void makePairList(std::vector< std::vector<size_t>>& pairs_list);

  void makePairToEmission(const std::vector< std::vector<size_t>>& pairs_list,
                          std::vector<std::vector<size_t> >& pairs_to_emission);

  void findEffectiveK();

  std::vector<size_t> makeEnumeratedArray();

  std::vector<size_t> activePairsArray(const std::vector<size_t>& pair_permute_row);

  static size_t nchoose2(size_t n);



};



}   // organization level namespace
}   // project level namespace



#endif
