
#include "kgd_panel.h"
#include <math.h>
#include <iostream>


namespace kgd = kellerberrin::deconvolv;


kgd::Panel::Panel() : TxtReader() {

  setTruePanelSize(0);
  setInbreedingPanelSize(0);

}


void kgd::Panel::readFromFile(const char inchar[]) {

  readFromFileBase(inchar);
  setTruePanelSize(getInfoLines());
  setInbreedingPanelSize(truePanelSize());

}


void kgd::Panel::checkForExceptions(size_t nLoci, std::string panelFileName) {

  if (getContent().size() != nLoci) {

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
  double nPanelMinus1 = nPanelDouble - 1.0;

  for (size_t i = 0; i < getPosition().size(); i++) {

    for (size_t j = 1; j < getPosition()[i].size(); j++) {

      geneticDistance = (getPositionIndexFloat(i,j) - getPositionIndexFloat(i, j - 1)) / averageMorganDistance;
      //rho = geneticDistance * 2 * Ne;
      rho = geneticDistance * G;

      double pRecTmp = (useConstRecomb) ? constRecombProb : 1.0 - exp(-rho);
      pRec_.push_back(pRecTmp);

      double pRecEachHapTmp = pRecTmp / nPanelDouble;
      pRecEachHap_.push_back(pRecTmp / nPanelDouble);

      double pNoRecTmp = 1.0 - pRecTmp;
      pNoRec_.push_back(pNoRecTmp);

      double secondPRecEachHapTmp = (forbidCopyFromSame) ? (pRecTmp / nPanelMinus1) : pRecEachHapTmp; // allowing copy from the same

      pRecRec_.push_back(pRecEachHapTmp * secondPRecEachHapTmp);
      pRecNoRec_.push_back(secondPRecEachHapTmp * pNoRecTmp);
      pNoRecNoRec_.push_back(pNoRecTmp * pNoRecTmp);

    }

    pRec_.push_back(1.0);
    pRecEachHap_.push_back(1.0 / nPanelDouble);
    pNoRec_.push_back(0.0);
    pRecRec_.push_back(
    ((forbidCopyFromSame) ? (1.0 / nPanelDouble / nPanelMinus1) : (1.0 / nPanelDouble / nPanelDouble)));
    pRecNoRec_.push_back(0.0);
    pNoRecNoRec_.push_back(0.0);

  }

  assert(pRec_.size() == getLoci());
  assert(pRecEachHap_.size() == getLoci());
  assert(pNoRec_.size() == getLoci());
  assert(pRecRec_.size() == getLoci());
  assert(pRecNoRec_.size() == getLoci());
  assert(pNoRecNoRec_.size() == getLoci());

}


void kgd::Panel::buildExamplePanel1() {

  addChrom("Pf3D7_01_v3");
  addPosition(std::vector<size_t>({93157, 94422, 94459, 94487, 95518, 95632, 95641}));
  addIndexOfChromStarts(0);
  buildExamplePanelContent();

}


void kgd::Panel::buildExamplePanel2() {

  addChrom("Pf3D7_01_v3");
  addChrom("Pf3D7_02_v3");
  addChrom("Pf3D7_03_v3");
  addPosition(std::vector<size_t>({93157}));
  addPosition(std::vector<size_t>({94422, 94459, 94487, 95518, 95632}));
  addPosition(std::vector<size_t>({95641}));
  addIndexOfChromStarts(0);
  addIndexOfChromStarts(1);
  addIndexOfChromStarts(6);
  buildExamplePanelContent();

}


void kgd::Panel::buildExamplePanelContent() {

  addContent(std::vector<double>({0, 0, 0, 1}));
  addContent(std::vector<double>({0, 0, 0, 1}));
  addContent(std::vector<double>({0, 0, 0, 1}));
  addContent(std::vector<double>({0, 0, 0, 1}));
  addContent(std::vector<double>({0, 1, 1, 0}));
  addContent(std::vector<double>({0, 0, 1, 0}));
  addContent(std::vector<double>({0, 0, 1, 0}));
  setLoci(getContent().size());
  setInfoLines(getContent().back().size());
  setTruePanelSize(getInfoLines());
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

  for (size_t siteI = 0; siteI < getContent().size(); siteI++) {

    for (size_t panelStrainJ = truePanelSize(); panelStrainJ < inbreedingPanelSize(); panelStrainJ++) {

      setContent()[siteI].push_back(1);

    }

    assert(inbreedingPanelSizeSetTo == getContent()[siteI].size());

  }

}


void kgd::Panel::updatePanelWithHaps(size_t inbreedingPanelSizeSetTo,
                                     size_t excludedStrain,
                                     const std::vector<std::vector<double> > &haps) {

  setInbreedingPanelSize(inbreedingPanelSizeSetTo);

  // If allows inbreeding, update reference panel_ by including strain haplotypes
  dout << "*************** Allow inbreeding **************" << std::endl;
  dout << "*** Update reference panel_ without strain " << excludedStrain << " ***" << std::endl;

  if (truePanelSize() == inbreedingPanelSize()) {

    return;

  }

  for (size_t siteI = 0; siteI < getContent().size(); siteI++) {

    size_t shiftAfter = inbreedingPanelSize();

    for (size_t panelStrainJ = truePanelSize(); panelStrainJ < inbreedingPanelSize(); panelStrainJ++) {

      size_t strainIndex = panelStrainJ - truePanelSize();

      if (strainIndex == excludedStrain) {

        shiftAfter = panelStrainJ;

      }

      if (shiftAfter <= panelStrainJ) {

        strainIndex++;

      }

      setContent()[siteI][panelStrainJ] = haps[siteI][strainIndex];

    }

  }

}

