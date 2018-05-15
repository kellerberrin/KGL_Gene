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

#include "kgd_panel.h"
#include <math.h>
#include <iostream>


namespace kgd = kellerberrin::deploid;


kgd::Panel::Panel() : TxtReader() {

  setTruePanelSize(0);
  setInbreedingPanelSize(0);

}


void kgd::Panel::readFromFile(const char inchar[]) {

  readFromFileBase(inchar);
  setTruePanelSize(nInfoLines_);
  setInbreedingPanelSize(truePanelSize());

}


void kgd::Panel::checkForExceptions(size_t nLoci, std::string panelFileName) {

  if (content_.size() != nLoci) {

    throw LociNumberUnequal(panelFileName);

  }

  if (pRec_.size() != nLoci) {

    throw LociNumberUnequal(panelFileName);

  }

  return;

}


void kgd::Panel::computeRecombProbs(double averageCentimorganDistance,
                                    double G,
                                    bool useConstRecomb,
                                    double constRecombProb,
                                    bool forbidCopyFromSame) {

  assert(pRec_.size() == 0);
  assert(pRecEachHap_.size() == 0);
  assert(pNoRec_.size() == 0);
  assert(pRecRec_.size() == 0);
  assert(pRecNoRec_.size() == 0);
  assert(pNoRecNoRec_.size() == 0);

  double averageMorganDistance = averageCentimorganDistance * 100;
  double geneticDistance;
  double rho;
  double nPanelDouble = (double) truePanelSize();
  double nPanlelMinus1 = nPanelDouble - 1.0;

  for (size_t i = 0; i < position_.size(); i++) {

    for (size_t j = 1; j < position_[i].size(); j++) {

      geneticDistance = (double) (position_[i][j] - position_[i][j - 1]) / averageMorganDistance;
      //rho = geneticDistance * 2 * Ne;
      rho = geneticDistance * G;

      double pRecTmp = (useConstRecomb) ? constRecombProb : 1.0 - exp(-rho);
      pRec_.push_back(pRecTmp);

      double pRecEachHapTmp = pRecTmp / nPanelDouble;
      pRecEachHap_.push_back(pRecTmp / nPanelDouble);

      double pNoRecTmp = 1.0 - pRecTmp;
      pNoRec_.push_back(pNoRecTmp);

      double secondPRecEachHapTmp = (forbidCopyFromSame) ? (pRecTmp / nPanlelMinus1)
                                                         : pRecEachHapTmp; // allowing copy from the same

      pRecRec_.push_back(pRecEachHapTmp * secondPRecEachHapTmp);
      pRecNoRec_.push_back(secondPRecEachHapTmp * pNoRecTmp);
      pNoRecNoRec_.push_back(pNoRecTmp * pNoRecTmp);

    }

    pRec_.push_back(1.0);
    pRecEachHap_.push_back(1.0 / nPanelDouble);
    pNoRec_.push_back(0.0);
    pRecRec_.push_back(
    ((forbidCopyFromSame) ? (1.0 / nPanelDouble / nPanlelMinus1) : (1.0 / nPanelDouble / nPanelDouble)));
    pRecNoRec_.push_back(0.0);
    pNoRecNoRec_.push_back(0.0);

  }

  assert(pRec_.size() == nLoci_);
  assert(pRecEachHap_.size() == nLoci_);
  assert(pNoRec_.size() == nLoci_);
  assert(pRecRec_.size() == nLoci_);
  assert(pRecNoRec_.size() == nLoci_);
  assert(pNoRecNoRec_.size() == nLoci_);

}


void kgd::Panel::buildExamplePanel1() {

  chrom_ = std::vector<std::string>({"Pf3D7_01_v3"});
  position_.push_back(std::vector<int>({93157, 94422, 94459, 94487, 95518, 95632, 95641}));
  indexOfChromStarts_ = std::vector<size_t>({0});
  buildExamplePanelContent();

}


void kgd::Panel::buildExamplePanel2() {

  chrom_ = std::vector<std::string>({"Pf3D7_01_v3", "Pf3D7_02_v3", "Pf3D7_03_v3"});
  position_.push_back(std::vector<int>({93157}));
  position_.push_back(std::vector<int>({94422, 94459, 94487, 95518, 95632}));
  position_.push_back(std::vector<int>({95641}));
  indexOfChromStarts_ = std::vector<size_t>({0, 1, 6});
  buildExamplePanelContent();

}


void kgd::Panel::buildExamplePanelContent() {

  content_.push_back(std::vector<double>({0, 0, 0, 1}));
  content_.push_back(std::vector<double>({0, 0, 0, 1}));
  content_.push_back(std::vector<double>({0, 0, 0, 1}));
  content_.push_back(std::vector<double>({0, 0, 0, 1}));
  content_.push_back(std::vector<double>({0, 1, 1, 0}));
  content_.push_back(std::vector<double>({0, 0, 1, 0}));
  content_.push_back(std::vector<double>({0, 0, 1, 0}));
  nLoci_ = content_.size();
  nInfoLines_ = content_.back().size();
  setTruePanelSize(nInfoLines_);
  setInbreedingPanelSize(truePanelSize());

}


void kgd::Panel::initializeUpdatePanel(size_t inbreedingPanelSizeSetTo) {

  setInbreedingPanelSize(inbreedingPanelSizeSetTo);

  // If allows inbreeding, update reference panel_ by including strain haplotypes
  dout << "************* Allow inbreeding ************" << std::endl;
  dout << "** Initialize inbreeding reference panel_ **" << std::endl;

  if (truePanelSize() == inbreedingPanelSize()) {

    return;

  }

  for (size_t siteI = 0; siteI < content_.size(); siteI++) {

    for (size_t panelStrainJ = truePanelSize(); panelStrainJ < inbreedingPanelSize(); panelStrainJ++) {

      content_[siteI].push_back(1);

    }

    assert(inbreedingPanelSizeSetTo == content_[siteI].size());

  }

}


void kgd::Panel::updatePanelWithHaps(size_t inbreedingPanelSizeSetTo,
                                     size_t excludedStrain,
                                     std::vector<std::vector<double> > &haps) {

  setInbreedingPanelSize(inbreedingPanelSizeSetTo);

  // If allows inbreeding, update reference panel_ by including strain haplotypes
  dout << "*************** Allow inbreeding **************" << std::endl;
  dout << "*** Update reference panel_ without strain " << excludedStrain << " ***" << std::endl;

  if (truePanelSize() == inbreedingPanelSize()) {

    return;

  }

  for (size_t siteI = 0; siteI < content_.size(); siteI++) {

    size_t shiftAfter = inbreedingPanelSize();

    for (size_t panelStrainJ = truePanelSize(); panelStrainJ < inbreedingPanelSize(); panelStrainJ++) {

      size_t strainIndex = panelStrainJ - truePanelSize();

      if (strainIndex == excludedStrain) {

        shiftAfter = panelStrainJ;

      }

      if (shiftAfter <= panelStrainJ) {

        strainIndex++;

      }

      content_[siteI][panelStrainJ] = haps[siteI][strainIndex];

    }

  }

}


void kgd::IBDrecombProbs::computeRecombProbs(double averageCentimorganDistance,
                                             double G,
                                             bool useConstRecomb,
                                             double constRecombProb) {

  assert(pRec_.size() == 0);
  assert(pNoRec_.size() == 0);

  double averageMorganDistance = averageCentimorganDistance * 100;
  double geneticDistance;
  double rho;

  for (size_t i = 0; i < position_.size(); i++) {

    for (size_t j = 1; j < position_[i].size(); j++) {

      geneticDistance = (double) (position_[i][j] - position_[i][j - 1]) / averageMorganDistance;
      //rho = geneticDistance * 2 * Ne;
      rho = geneticDistance * G;
      double pRecTmp = (useConstRecomb) ? constRecombProb : 1.0 - exp(-rho);
      pRec_.push_back(pRecTmp);
      double pNoRecTmp = 1.0 - pRecTmp;
      pNoRec_.push_back(pNoRecTmp);

    }

    pRec_.push_back(1.0);
    pNoRec_.push_back(0.0);

  }

  assert(pRec_.size() == nLoci_);
  assert(pNoRec_.size() == nLoci_);

}



//vector<vector<double>> outtrans(out[0].size(),
//vector<double>(out.size()));
//for (size_t i = 0; i < out.size(); ++i)
//for (size_t j = 0; j < out[0].size(); ++j)
//outtrans[j][i] = out[i][j];
