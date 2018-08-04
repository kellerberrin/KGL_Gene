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


#ifndef KGD_HAP_H
#define KGD_HAP_H


#include <vector>
#include <iostream>
#include "kgd_utility.h"
#include "kgd_panel.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class UpdateHap {
#ifdef UNITTEST
  friend class TestUpdatePairHap;
  friend class TestUpdateSingleHap;
  friend class TestUpdateHap;
#endif


public:

  UpdateHap(size_t kStrain,
            size_t segmentStartIndex,
            size_t nLoci,
            std::shared_ptr<Panel> panel,
            double missCopyProb,
            double scalingFactor);

  virtual ~UpdateHap() = default;

  // Access functions.
  size_t nPanel() const { return nPanel_; }
  double getNewLLKIndex(size_t index) const { return newLLK[index]; }

  // Modification functions.
  void setPanelSize(const size_t setTo) { nPanel_ = setTo; }

protected:

  std::shared_ptr<Panel> panel_;
  double missCopyProb_;
  size_t kStrain_;
  size_t nPanel_;
  std::vector<double> newLLK;
  size_t segmentStartIndex_;
  size_t nLoci_;
  std::vector<std::vector<double> > emission_;
  double scalingFactor_;

  double scalingFactor() const { return scalingFactor_; }
  void setScalingFactor(const double setTo) { scalingFactor_ = setTo; }


};




}   // organization level namespace
}   // project level namespace



#endif
