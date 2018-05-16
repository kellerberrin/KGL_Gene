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

#include "kgd_deconvolv_app.h"
#include "kgd_exceptions.h"
#include "kgd_txt_reader.h"
#include <algorithm> // find
#include <iostream>


namespace kgd = kellerberrin::deconvolv;


kgd::VariantIndex::VariantIndex() {

  init();

}

void kgd::VariantIndex::findWhoToBeKept(std::shared_ptr<ExcludeMarker> excludedMarkers) {

  dout << " Starts findWhoToBeKept " << std::endl;

  assert (indexOfContentToBeKept.size() == 0);
  assert (indexOfPosToBeKept.size() == 0);

  for (size_t chromI = 0; chromI < chrom_.size(); chromI++) {

    dout << "   Going through chrom " << chrom_[chromI] << std::endl;

    std::vector<size_t> tmpindexOfPosToBeKept;

    // detemine if something needs to be removed from the current chrom.
    auto chromIt = find(excludedMarkers->chrom_.begin(), excludedMarkers->chrom_.end(), chrom_[chromI]);

    size_t hapIndex = indexOfChromStarts_[chromI];

    size_t chromIndexInExclude = std::distance(excludedMarkers->chrom_.begin(), chromIt);

    for (size_t posI = 0; posI < position_[chromI].size(); posI++) {

      if (chromIt == excludedMarkers->chrom_.end()) {

        indexOfContentToBeKept.push_back(hapIndex);
        tmpindexOfPosToBeKept.push_back(posI);

      } else if (std::find(excludedMarkers->position_[chromIndexInExclude].begin(),
                           excludedMarkers->position_[chromIndexInExclude].end(),
                           position_[chromI][posI]) == excludedMarkers->position_[chromIndexInExclude].end()) {

        indexOfContentToBeKept.push_back(hapIndex);
        tmpindexOfPosToBeKept.push_back(posI);

      }

      hapIndex++;

    }

    indexOfPosToBeKept.push_back(tmpindexOfPosToBeKept);

  }

  assert (indexOfPosToBeKept.size() == chrom_.size());

  dout << indexOfContentToBeKept.size() << " sites need to be Kept " << std::endl;

}


void kgd::VariantIndex::findAndKeepMarkers(std::shared_ptr<ExcludeMarker> excludedMarkers) {

  setDoneGetIndexOfChromStarts(false);

  dout << " findAndKeepMarkers called" << std::endl;

  findWhoToBeKept(excludedMarkers);
  removePositions();
  IndexOfChromStarts();
  removeMarkers();

}


void kgd::VariantIndex::removePositions() {

  assert (keptPosition_.size() == (size_t) 0);

  for (size_t chromI = 0; chromI < chrom_.size(); chromI++) {

    std::vector<int> tmpKeptPosition_;

    for (size_t i = 0; i < indexOfPosToBeKept[chromI].size(); i++) {

      tmpKeptPosition_.push_back(position_[chromI][indexOfPosToBeKept[chromI][i]]);

    }

    keptPosition_.push_back(tmpKeptPosition_);

  }

  position_.clear();
  position_ = keptPosition_;
  keptPosition_.clear();

}


void kgd::VariantIndex::IndexOfChromStarts() {

  assert(not doneGetIndexOfChromStarts_);

  indexOfChromStarts_.clear();

  assert(indexOfChromStarts_.size() == 0);

  indexOfChromStarts_.push_back((size_t) 0);

  for (size_t tmpChrom = 0; indexOfChromStarts_.size() < chrom_.size(); tmpChrom++) {

    size_t chrom_start = indexOfChromStarts_.back() + position_[tmpChrom].size();

    indexOfChromStarts_.push_back(chrom_start);

    ExecEnv::log().info("IndexOfChromStarts(); Chrom: {}, positions: {}, chrom_start: {}", tmpChrom, position_[tmpChrom].size(), chrom_start);

  }

  assert(indexOfChromStarts_.size() == chrom_.size());



  setDoneGetIndexOfChromStarts(true);

}


void kgd::VariantIndex::init() {

  setDoneGetIndexOfChromStarts(false);

}

void kgd::VariantIndex::removeMarkers() { throw VirtualFunctionShouldNotBeCalled(); }
