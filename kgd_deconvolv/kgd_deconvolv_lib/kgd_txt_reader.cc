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

#include <fstream>
#include <iostream>
#include <iterator>     // std::distance
#include "kgd_exceptions.h"
#include "kgd_txt_reader.h"


namespace kgd = kellerberrin::deconvolv;


void kgd::TxtReader::readFromFileBase(const char inchar[]) {

  fileName_ = std::string(inchar);
  tmpChromInex_ = -1;

  std::ifstream in_file(inchar);
  std::string tmp_line;

  if (in_file.good()) {

    getline(in_file, tmp_line); // skip the first line, which is the header
    getline(in_file, tmp_line);

    while (tmp_line.size() > 0) {

      size_t field_start = 0;
      size_t field_end = 0;
      size_t field_index = 0;
      std::vector<double> contentRow;

      while (field_end < tmp_line.size()) {

        field_end = std::min(std::min(tmp_line.find(',', field_start),
                                      tmp_line.find('\t', field_start)),
                             tmp_line.find('\n', field_start));

        std::string tmp_str = tmp_line.substr(field_start, field_end - field_start);

        if (field_index > 1) {

          contentRow.push_back(strtod(tmp_str.c_str(), NULL));

        } else if (field_index == 0) {

          extractChrom(tmp_str);

        } else if (field_index == 1) {

          extractPOS(tmp_str);

        }

        field_start = field_end + 1;
        field_index++;

      }

      content_.push_back(contentRow);
      getline(in_file, tmp_line);

    }

  } else {

    throw InvalidInputFile(fileName_);

  }

  in_file.close();

  allele_reader_.parseFile(fileName_);

  addPosition(tmpPosition_);

  setLoci(content_.size());
  nInfoLines_ = content_.back().size();

  if (nInfoLines_ == 1) {

    reshapeContentToInfo();

  }

  IndexOfChromStarts();

  assert (tmpChromInex_ > -1);
  assert (getChrom().size() == getPosition().size());
  assert(getDoneIndexOfChromStarts());

}


void kgd::TxtReader::extractChrom(std::string &tmp_str) {

  if (tmpChromInex_ >= 0) {

    if (tmp_str != getChrom().back()) {

      tmpChromInex_++;
      // save current positions
      addPosition(tmpPosition_);

      // start new chrom
      tmpPosition_.clear();
      addChrom(tmp_str);

    }

  } else {

    tmpChromInex_++;

    assert (getChrom().size() == 0);

    addChrom(tmp_str);

    assert (tmpPosition_.size() == 0);
    assert (getPosition().size() == 0);

  }

}


void kgd::TxtReader::extractPOS(std::string &tmp_str) {
  int ret;

  try {

    ret = std::stoi(tmp_str.c_str(), NULL);

  } catch (const std::exception &e) {

    throw BadConversion(tmp_str, fileName_);
  }

  tmpPosition_.push_back(ret);

}


void kgd::TxtReader::reshapeContentToInfo() {

  assert (info_.size() == 0);

  for (size_t i = 0; i < content_.size(); i++) {

    info_.push_back(content_[i][0]);

  }

}


void kgd::TxtReader::removeMarkers() {

  assert(keptContent_.size() == 0);

  for (auto const &value: getIndexOfContentToBeKept()) {

    keptContent_.push_back(content_[value]);

  }

  content_.clear();
  assert(content_.size() == (size_t) 0);
  content_ = keptContent_;
  keptContent_.clear();

  if (nInfoLines_ == 1) {

    info_.clear();
    reshapeContentToInfo();

  }

  setLoci(content_.size());

}
