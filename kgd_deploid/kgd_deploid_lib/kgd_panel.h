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

#ifndef KGD_PANEL_H
#define KGD_PANEL_H


#include "kgd_txt_reader.h"
#include "kgd_exceptions.h"



namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



class Panel : public TxtReader {
#ifdef UNITTEST
  friend class TestPanel;
  friend class TestInitialHaplotypes;
  friend class TestUpdateHap;
  friend class TestUpdatePairHap;
  friend class TestUpdateSingleHap;
#endif

  friend class McmcMachinery;

  friend class UpdateSingleHap;

  friend class UpdatePairHap;

  friend class UpdateHap;

  friend class DEploidIO;

  friend class InitialHaplotypes;

public:

  Panel();
  virtual ~Panel() {};

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


  void setTruePanelSize(const size_t setTo) { this->truePanelSize_ = setTo; }

  void setInbreedingPanelSize(const size_t setTo) { this->inbreedingPanelSize_ = setTo; }

  size_t inbreedingPanelSize() const { return this->inbreedingPanelSize_; }

  size_t truePanelSize() const { return this->truePanelSize_; }

  // Methods
  void readFromFile(const char inchar[]);

  void computeRecombProbs(double averageCentimorganDistance,
                          double Ne,
                          bool useConstRecomb,
                          double constRecombProb,
                          bool forbidCopyFromSame);

  void checkForExceptions(size_t nLoci, std::string panelFileName);

  void initializeUpdatePanel(size_t inbreedingPanelSizeSetTo);

  void updatePanelWithHaps(size_t inbreedingPanelSizeSetTo,
                           size_t excludedStrain,
                           std::vector<std::vector<double> > &haps);

  void print();

  void buildExamplePanelContent();

  void buildExamplePanel1();

  void buildExamplePanel2();

};


class InitialHaplotypes : public Panel {

  #ifdef UNITTEST
  friend class TestInitialHaplotypes;
#endif

  friend class DEploidIO;

  InitialHaplotypes() : Panel() {}
  ~InitialHaplotypes() = default;

};


class IBDrecombProbs : public VariantIndex {

  friend class IBDpath;

#ifdef UNITTEST
  friend class TestIBDpath;
#endif

private:

  std::vector<double> pRec_;
  std::vector<double> pNoRec_; // = 1.0 - pRec;

  void computeRecombProbs(double averageCentimorganDistance, double Ne, bool useConstRecomb, double constRecombProb);

public:

  IBDrecombProbs() = default;
  IBDrecombProbs(std::vector<std::vector<int> > position, size_t nLoci) {

    this->position_ = position;
    this->nLoci_ = nLoci;

  }

  ~IBDrecombProbs() = default;

};



}   // organization level namespace
}   // project level namespace



#endif
