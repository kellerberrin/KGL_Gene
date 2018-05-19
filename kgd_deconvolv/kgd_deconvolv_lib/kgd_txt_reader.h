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


#ifndef KGD_TXTREADER_H
#define KGD_TXTREADER_H


#include "kgd_variant_index.h"
#include "kgd_exceptions.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


class TxtReader : public VariantIndex {

#ifdef UNITTEST
  friend class TestPanel;
  friend class TestTxtReader;
  friend class TestInitialHaplotypes;
#endif


public: // move the following to private

  TxtReader() = default;
  ~TxtReader() override = default;

  // Methods
  virtual void readFromFile(const char inchar[]) { this->readFromFileBase(inchar); };
  void readFromFileBase(const char inchar[]);
  void removeMarkers();

  // Access routines.
  double getContentIndex(size_t loci, size_t strain) const { return content_[loci][strain]; }
  const std::vector<std::vector<double> >& getContent() const { return content_; }
  size_t getInfoLines() const { return nInfoLines_; }
  const std::vector<double>& getInfo() const { return info_; }

  void setInfoLines(size_t lines) { nInfoLines_ = lines; }
  void addContent(const std::vector<double>& content) { content_.push_back(content); }
  std::vector<std::vector<double> >& setContent() { return content_; }

private:

  // content is a matrix of n.loci by n.strains, i.e. content length is n.loci
  std::vector<std::vector<double> > content_;
  std::vector<std::vector<double> > keptContent_;
  // info_ only refers to the first column of the content
  std::vector<double> info_;

  std::vector<int> tmpPosition_;

  size_t nInfoLines_;
  int tmpChromInex_;

  std::string fileName_;

  // Methods
  void extractChrom(std::string &tmp_str);
  void extractPOS(std::string &tmp_str);
  void reshapeContentToInfo();


};


class ExcludeMarker : public TxtReader {

  // sorting

public:

  ExcludeMarker() = default;
  ~ExcludeMarker() override = default;

};


}   // organization level namespace
}   // project level namespace


#endif
