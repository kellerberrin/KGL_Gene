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

#ifndef KGD_VARIANTINDEX_H
#define KGD_VARIANTINDEX_H


#include <vector>
#include <memory>
#include <string>
#include <cassert>
#include "kgd_global.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



// Forward Def.
class ExcludeMarker;


class VariantIndex {

public:

  VariantIndex();
  virtual ~VariantIndex() = default;

  // Access functions.
  const std::vector<std::string>& getChrom() const { return chrom_; }
  const std::vector<size_t>& getIndexChromStarts() const { return indexOfChromStarts_; }
  const std::vector<std::vector<int> >& getPosition() const { return position_; }
  double getPositionIndexFloat(size_t chrom, size_t idx) const { return static_cast<double>(position_[chrom][idx]); }
  bool getDoneIndexOfChromStarts() const { return doneGetIndexOfChromStarts_; }
  const std::vector<size_t>& getIndexOfContentToBeKept() const { return indexOfContentToBeKept; }
  size_t getLoci() const { return nLoci_; }

  void setLoci(size_t loci) { nLoci_ = loci; }
  void setPosition(const std::vector<std::vector<int>>& position) { position_ = position; }
  void addPosition(const std::vector<int>& position) { position_.push_back(position); }
  void addChrom(const std::string& chrom) { chrom_.push_back(chrom); }
  void addIndexOfChromStarts(size_t start) { indexOfChromStarts_.push_back(start); }

  // Methods
  void IndexOfChromStarts();
  void findAndKeepMarkers(std::shared_ptr<ExcludeMarker> excludedMarkers);
  void InitChrom() { chrom_.clear(); }
  void InitPosition() { position_.clear(); }

#ifdef UNITTEST
  friend class TestPanel;
  friend class TestTxtReader;
  friend class TestInitialHaplotypes;
#endif

private:

  std::vector<std::string> chrom_;
  std::vector<size_t> indexOfChromStarts_;
  std::vector<std::vector<int> > position_;
  std::vector<std::vector<int> > keptPosition_;
  size_t nLoci_;

  /* Index of content/info will be kept */
  std::vector<size_t> indexOfContentToBeKept;
  /* Index of positions entry to be kept, this will have the same size as this->chrom_, */
  std::vector<std::vector<size_t> > indexOfPosToBeKept;

  bool doneGetIndexOfChromStarts_;

  // Getters and Setters.
  void setDoneGetIndexOfChromStarts(bool setTo) { this->doneGetIndexOfChromStarts_ = setTo; }

  // Methods
  void init();
  void removePositions();
  // For removing markers and positions
  void findWhoToBeKept(std::shared_ptr<ExcludeMarker> excludedMarkers);

  virtual void removeMarkers();


};


}   // organization level namespace
}   // project level namespace



#endif
