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

#ifndef KGD_PANEL_H
#define KGD_PANEL_H


#include "kgd_txt_reader.h"
#include "kgd_exceptions.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class Panel : public TxtReader {

#ifdef UNITTEST
  friend class TestPanel;
  friend class TestInitialHaplotypes;
  friend class TestUpdateHap;
  friend class TestUpdatePairHap;
  friend class TestUpdateSingleHap;
#endif

public:

  Panel();
  ~Panel() override = default;

  // Access functions.
  size_t truePanelSize() const { return this->truePanelSize_; }
  size_t inbreedingPanelSize() const { return this->inbreedingPanelSize_; }
  const std::vector<double>& getRec() const { return pRec_; }
  double getRecIndex(size_t index) const { return pRec_[index]; }
  double getRecEachHapIndex(size_t index) const { return pRecEachHap_[index]; }
  double getNoRecIndex(size_t index) const { return pNoRec_[index]; }
  double getRecRecIndex(size_t index) const { return pRecRec_[index]; }
  double getRecNoRecIndex(size_t index) const { return pRecNoRec_[index]; }
  double getNoRecNoRecIndex(size_t index) const { return pNoRecNoRec_[index]; }

  // Methods
  void initializeUpdatePanel(size_t inbreedingPanelSizeSetTo);
  void updatePanelWithHaps(size_t inbreedingPanelSizeSetTo,
                           size_t excludedStrain,
                           const std::vector<std::vector<double> > &haps);
  void readFromFile(const char inchar[]);
  void checkForExceptions(size_t nLoci, std::string panelFileName);
  void computeRecombProbs(double averageCentimorganDistance,
                          double Ne,
                          bool useConstRecomb,
                          double constRecombProb,
                          bool forbidCopyFromSame);

private:

  // Members
  std::vector<double> pRec_;
  // Used in update single haplotype
  std::vector<double> pRecEachHap_; // = pRec / nPanel_;
  std::vector<double> pNoRec_; // = 1.0 - pRec;
  // Used in update pair of haplotypes
  std::vector<double> pRecRec_; // pRecEachHap * pRecEachHap;
  std::vector<double> pRecNoRec_; // pRecEachHap * pNoRec;
  std::vector<double> pNoRecNoRec_; // pNoRec * pNoRec;

  size_t truePanelSize_;
  size_t inbreedingPanelSize_;

  // Setters and Getters.
  void setTruePanelSize(const size_t setTo) { this->truePanelSize_ = setTo; }
  void setInbreedingPanelSize(const size_t setTo) { this->inbreedingPanelSize_ = setTo; }

  // Methods
  void print();
  void buildExamplePanelContent();
  void buildExamplePanel1();
  void buildExamplePanel2();

};



}   // organization level namespace
}   // project level namespace


#endif
