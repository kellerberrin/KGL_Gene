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
*/

#include <ctime>
#include <iterator>
#include <cassert>        // assert
#include <iomanip>        // std::setw
#include "kgd_utility.h"    // normailize by sum
#include "kgd_update_haplotype.h"  // chromPainting
#include "kgd_update_single_haplotype.h"  // chromPainting
#include "kgd_deploid_io.h"
#include "kgd_ibdpath.h"
#include "kgd_deconvolv_app.h"

namespace kgd = kellerberrin::deconvolv;


kgd::DEploidIO::DEploidIO() {

  init(); // Reset to default values before parsing
  parse();
  checkInput();
  readFiles();
  finalize();
  removeFilesWithSameName();

}


kgd::DEploidIO::DEploidIO(const MixtureDataObj& mixture_data) {

  init(); // Reset to default values before parsing
  setMixtureData(mixture_data);
  finalize();
  removeFilesWithSameName();

}


void kgd::DEploidIO::init() {

  getTime(true);

  panel_ = nullptr;
  initialProp_.clear();

  setkStrain(5);

  refFileName_.clear();
  altFileName_.clear();
  plafFileName_.clear();
  panelFileName_.clear();
  excludeFileName_.clear();

  prefix_ = "kgd_deconvolv";
  compileTime_ = "";
  dEploidGitVersion_ = Deconvolv::VERSION;

}





void kgd::DEploidIO::getTime(bool isStartingTime) {

  time_t now = time(0);
  // convert now to string form

  char *dt = ctime(&now);

  if (isStartingTime) {

    startingTime_ = dt;

  } else {

    endTime_ = dt;

  }

}

void kgd::DEploidIO::readFiles() {

  if (getMixtureControl().useVcf()) { // read vcf files, and parse it to refCount and altCount

    if (getMixtureControl().excludeSites()) {

      mixture_data_.readVCFPlafExclude(vcfFileName_, plafFileName_, excludeFileName_);

    } else {

      mixture_data_.readVCFPlaf(vcfFileName_, plafFileName_);

    }

  } else {

    if (getMixtureControl().excludeSites()) {

      mixture_data_.readRefAltPlafExclude(refFileName_, altFileName_, plafFileName_, excludeFileName_);

    } else {

      mixture_data_.readRefAltPlaf(refFileName_, altFileName_, plafFileName_);

    }

  }

  readPanel();

}

void kgd::DEploidIO::finalize() {

  if (getMixtureControl().doIbdPainting() or getMixtureControl().doComputeLLK()) {

    if (!initialPropWasGiven()) {

      throw InitialPropUngiven("");

    }

  }

  if (useIBD() && kStrain() == 1) {

    throw InvalidK();

  }

  if (getMixtureControl().compressVcf() && !getMixtureControl().doExportVcf()) {

    throw VcfOutUnSpecified("");

  }

  removeFilesWithSameName();

  IBDpathChangeAt_ = std::vector<double>(nLoci());
  finalIBDpathChangeAt_ = std::vector<double>(nLoci());

  siteOfTwoSwitchOne_ = std::vector<double>(nLoci());
  siteOfTwoMissCopyOne_ = std::vector<double>(nLoci());
  siteOfTwoSwitchTwo_ = std::vector<double>(nLoci());
  siteOfTwoMissCopyTwo_ = std::vector<double>(nLoci());
  siteOfOneSwitchOne_ = std::vector<double>(nLoci());
  siteOfOneMissCopyOne_ = std::vector<double>(nLoci());

  finalSiteOfTwoSwitchOne_ = std::vector<double>(nLoci());
  finalSiteOfTwoMissCopyOne_ = std::vector<double>(nLoci());
  finalSiteOfTwoSwitchTwo_ = std::vector<double>(nLoci());
  finalSiteOfTwoMissCopyTwo_ = std::vector<double>(nLoci());
  finalSiteOfOneSwitchOne_ = std::vector<double>(nLoci());
  finalSiteOfOneMissCopyOne_ = std::vector<double>(nLoci());

}


void kgd::DEploidIO::removeFilesWithSameName() {

  strExportProp_ = prefix_ + ".prop";
  strExportLLK_ = prefix_ + ".llk";
  strExportHap_ = prefix_ + ".hap";

  strIbdExportProp_ = prefix_ + ".ibd.prop";
  strIbdExportLLK_ = prefix_ + ".ibd.llk";
  strIbdExportHap_ = prefix_ + ".ibd.hap";
  strIbdExportProbs_ = prefix_ + ".ibd.probs";

  strExportVcf_ = prefix_ + ".vcf";

  if (getMixtureControl().compressVcf()) {

    strExportVcf_ += ".gz";

  }

  strExportLog_ = prefix_ + ((getMixtureControl().doLsPainting()) ? ".painting" : "") + ".log";
  strExportRecombProb_ = prefix_ + ".recomb";

  strExportExtra_ = prefix_ + ".extra";

  if (getMixtureControl().doLsPainting()) {

    if (useIBD()) {

      remove(strIbdExportProp_.c_str());
      remove(strIbdExportLLK_.c_str());
      remove(strIbdExportHap_.c_str());

    }

    remove(strExportLLK_.c_str());
    remove(strExportHap_.c_str());
    remove(strExportVcf_.c_str());
    remove(strExportProp_.c_str());
    remove(strExportExtra_.c_str());
    remove(strIbdExportProbs_.c_str());

  }

  if (getMixtureControl().doLsPainting() || getMixtureControl().doExportPostProb()) {

    if (useIBD()) {

      strIbdExportSingleFwdProbPrefix_ = prefix_ + ".ibd.single";

      for (size_t i = 0; i < kStrain(); i++) {

        std::string tmpStrExportSingleFwdProb = strIbdExportSingleFwdProbPrefix_ + std::to_string(i);
        remove(tmpStrExportSingleFwdProb.c_str());

      }

      strIbdExportPairFwdProb_ = prefix_ + ".ibd.pair";
      remove(strIbdExportPairFwdProb_.c_str());

    }

    strExportSingleFwdProbPrefix_ = prefix_ + ".single";

    for (size_t i = 0; i < kStrain(); i++) {

      std::string tmpStrExportSingleFwdProb = strExportSingleFwdProbPrefix_ + std::to_string(i);
      remove(tmpStrExportSingleFwdProb.c_str());

    }

    strExportPairFwdProb_ = prefix_ + ".pair";
    remove(strExportPairFwdProb_.c_str());
  }

  remove(strExportLog_.c_str());
  remove(strExportRecombProb_.c_str());

}


void kgd::DEploidIO::parse() {

  const kgd::DeconvolvArgs &args = Deconvolv::getArgs();

  if (args.vcfFile!= kgd::DeconvolvArgs::NOT_SPECIFIED) {  // -vcfOut

    vcfFileName_ = args.vcfFile; // -vcf
    getMixtureControl().setUseVcf(true);

  } else {

    if (args.refFile!= kgd::DeconvolvArgs::NOT_SPECIFIED && args.altFile!= kgd::DeconvolvArgs::NOT_SPECIFIED) {

      refFileName_  = args.refFile;
      altFileName_ = args.altFile;
      getMixtureControl().setUseVcf(false);

    } else {

      ExecEnv::log().warn("No VCF '-vcf' or ALT and REF files '-ref', '-alt' specified. Assuming direct call.");

    }

  }

  if (args.vcfOutFile != kgd::DeconvolvArgs::NOT_SPECIFIED) {  // -vcfOut

    getMixtureControl().setDoExportVcf(true);

  }

  plafFileName_ = args.plafFile;  // -plaf

  if (args.panelFile != kgd::DeconvolvArgs::NOT_SPECIFIED) { // -panel_

    panelFileName_ = args.panelFile;
    getMixtureControl().setUsePanel(true);

  } else {

    getMixtureControl().setUsePanel(false);

  }

  getMixtureControl().setUsePanel(not args.noPanelFlag);

  if (not usePanel()) {

    getMixtureControl().setDoExportPostProb(false);
    getMixtureControl().setDoExportSwitchMissCopy(false);
    getMixtureControl().setDoAllowInbreeding(false);

  }

  if (args.excludeFile != kgd::DeconvolvArgs::NOT_SPECIFIED) {  // -exclude

    getMixtureControl().setExcludeSites(true);
    excludeFileName_ = args.excludeFile;

  }

  prefix_ = args.outputTemplate; // -o
  setkStrain(args.maxStrains); // -k
  getMixtureControl().setKStrainWasManuallySet(true);
  hapParameters().setMcmcSample(args.MCMCSamples); // -nSamples
  hapParameters().setMcmcBurn(args.MCMCBurnRate);  // -burn
  hapParameters().setMcmcMachineryRate(args.MCMCSampleRate);  // -rate

  if (args.inbreedingProbabilitiesFlag) { // -inbreeding

    getMixtureControl().setDoAllowInbreeding(true);
    getMixtureControl().setDoExportPostProb(true);

  }

  if (args.identityByDescentPainting) getMixtureControl().setDoIbdPainting(true); // -idbPainting

  setInitialProp(args.initialStrainProportions); // -initialP
  getMixtureControl().setUseIBD(args.identityByDescentFlag); // -idb

}

void kgd::DEploidIO::checkInput() {

  if (refFileName_.empty() && not getMixtureControl().useVcf()) {

    throw FileNameMissing("Ref count");

  }
  if (altFileName_.empty() && not getMixtureControl().useVcf()) {

    throw FileNameMissing("Alt count");

  }
  if (plafFileName_.empty()) {

    throw FileNameMissing("PLAF");
  }
  if (usePanel() && panelFileName_.empty() && not getMixtureControl().doIbdPainting() && not getMixtureControl().doComputeLLK()) {

    throw FileNameMissing("Reference panel_");

  }
  if (initialPropWasGiven() && (abs(Utility::sumOfVec(initialProp_) - 1.0) > 0.00001)) {

    throw SumOfPropNotOne(std::to_string(Utility::sumOfVec(initialProp_)));

  }

}


std::vector<double> kgd::DEploidIO::computeExpectedWsafFromInitialHap() {
  // Make this a separate function
  // calculate expected wsaf

  std::vector<double> expectedWsaf(getInitialHap().size(), 0.0);

  for (size_t i = 0; i < getInitialHap().size(); i++) {

    assert(kStrain() == getInitialHap()[i].size());

    for (size_t k = 0; k < kStrain(); k++) {

      expectedWsaf[i] += getInitialHap()[i][k] * finalProp_[k];

    }

    assert (expectedWsaf[i] >= 0);

  }

  return expectedWsaf;

}


void kgd::DEploidIO::computeLLKfromInitialHap() {

  for (auto const &value: initialProp_) {

    finalProp_.push_back(value);

  }

  std::vector<double> expectedWsaf = computeExpectedWsafFromInitialHap();

  if (expectedWsaf.size() != getMixtureData().getRefCount().size()) {

    throw LociNumberUnequal("Hap length differs from data!");

  }

  std::vector<double> llk = Utility::calcLLKs(getMixtureData().getRefCount(),
                                              getMixtureData().getAltCount(),
                                              expectedWsaf,
                                              0,
                                              expectedWsaf.size(),
                                              hapParameters().proposalUpdateScaling());

  llkFromInitialHap_ = Utility::sumOfVec(llk);

}


void kgd::DEploidIO::chromPainting() {

  ExecEnv::log().info("Painting haplotypes in :{}", initialHapFileName_);

  if (initialPropWasGiven() == false) {

    ExecEnv::log().info("Initial proportion was not specified. Set even proportions");

    double evenProp = 1.0 / (double) kStrain();

    for (size_t i = 0; i < kStrain(); i++) {

      initialProp_.push_back(evenProp);

    }

  }

  for (auto const &value: initialProp_) {

    finalProp_.push_back(value);

  }

  /// Painting posterior probabilities
  /// Export the p'
  /// Make this a separate class
  /// vector < vector <double> > hap;
  /// for ( size_t siteI = 0; siteI < decovolutedStrainsToBeRead.content_.size(); siteI++ ){
  /// vector <double> tmpHap;
  /// for ( size_t tmpk = 0; tmpk < kStrain_; tmpk++ ){
  /// tmpHap.push_back(decovolutedStrainsToBeRead.content_[siteI][tmpk]);
  /// }
  /// hap.push_back(tmpHap);
  /// }

  //vector < vector <double>> hap = decovolutedStrainsToBeRead.content_;

  std::vector<double> expectedWsaf = computeExpectedWsafFromInitialHap();

  if (doAllowInbreeding() == true) {

    panel_->initializeUpdatePanel(panel_->truePanelSize() + kStrain() - 1);

  }

  for (size_t tmpk = 0; tmpk < kStrain(); tmpk++) {

    if (doAllowInbreeding()) {

      panel_->updatePanelWithHaps(panel_->truePanelSize() + kStrain() - 1, tmpk, getInitialHap());

    }

    for (size_t chromi = 0; chromi < getMixtureData().indexOfChromStarts().size(); chromi++) {

      size_t start = getMixtureData().indexOfChromStarts()[chromi];
      size_t length = getMixtureData().getPosition()[chromi].size();

      ExecEnv::log().info("Painting Chrom: {} from site: {} to: {}", chromi, start, start + length);

      UpdateSingleHap updatingSingle(start,
                                     length,
                                     kStrain(),
                                     panel_,
                                     hapParameters().getMissCopyProb(),
                                     hapParameters().proposalUpdateScaling(),
                                     tmpk);

      if (doAllowInbreeding()) {

        updatingSingle.setPanelSize(panel_->inbreedingPanelSize());

      }

      updatingSingle.painting(getMixtureData().getRefCount(), getMixtureData().getAltCount(), expectedWsaf, finalProp_, initialHap_);

      writeLastSingleFwdProb(updatingSingle.getFwdBwdProbs(), chromi, tmpk, false); // false as not using ibd

    }

  }

}


void kgd::DEploidIO::readPanel() {

  if (usePanel() == false) {

    return;

  }

  if (getMixtureControl().doIbdPainting() or getMixtureControl().doComputeLLK()) {

    return;

  }

  panel_ = std::make_shared<Panel>();
  panel_->readFromFile(panelFileName_.c_str());

  if (getMixtureControl().excludeSites()) {

    std::shared_ptr<ExcludeMarker> excluded_reader_ptr(std::make_shared<ExcludeMarker>());
    excluded_reader_ptr->readFromFile(excludeFileName_.c_str());
    panel_->findAndKeepMarkers(excluded_reader_ptr);

  }

  panel_->computeRecombProbs(hapParameters().averageCentimorganDistance(),
                             hapParameters().parameterG(),
                             useConstRecomb(),
                             hapParameters().constRecombProb(),
                             forbidCopyFromSame());

  panel_->checkForExceptions(nLoci(), panelFileName_);

}


void kgd::DEploidIO::getIBDprobsIntegrated(const std::vector<std::vector<double> > &prob) {

  if (prob.size() != nLoci()) {

    throw InvalidInput("Invalid probabilities! Check size!");

  }

  assert(ibdProbsIntegrated_.size() == 0);

  for (size_t i = 0; i < prob[0].size(); i++) {

    ibdProbsIntegrated_.push_back(0.0);

  }

  for (size_t siteIndex = 0; siteIndex < nLoci(); siteIndex++) {

    for (size_t i = 0; i < prob[siteIndex].size(); i++) {

      ibdProbsIntegrated_[i] += prob[siteIndex][i];

    }

  }

  Utility::normalizeBySum(ibdProbsIntegrated_);

}


void kgd::DEploidIO::computeEffectiveKstrain(std::vector<double> proportion) {

  double tmpSumSq = 0.0;

  for (double p : proportion) {

    tmpSumSq += p * p;

  }

  effectiveKstrain_ = 1.0 / tmpSumSq;

}


void kgd::DEploidIO::computeInferredKstrain(std::vector<double> proportion) {

  inferredKstrain_ = 0;

  for (double p : proportion) {

    if (p > 0.01) {

      inferredKstrain_ += 1;

    }

  }

}


void kgd::DEploidIO::computeAdjustedEffectiveKstrain() {

  adjustedEffectiveKstrain_ = effectiveKstrain_;

  if ((inferredKstrain_ == 2) & (ibdProbsIntegrated_.size() == 2)) {

    if (ibdProbsIntegrated_[1] > 0.95) {

      adjustedEffectiveKstrain_ = 1;

    }

  }

}


void kgd::DEploidIO::paintIBD() {

  std::vector<double> goodProp;
  std::vector<size_t> goodStrainIdx;

  if (getMixtureControl().doIbdPainting()) {

    finalProp_ = initialProp_;

  }

  for (size_t i = 0; i < finalProp_.size(); i++) {

    if (finalProp_[i] > 0.01) {
      goodProp.push_back(finalProp_[i]);
      goodStrainIdx.push_back(i);
    }

  }

  if (goodProp.size() == 1) {

    return;

  }

  DEploidIO tmpDEploidIO(getMixtureData()); // (*this);
  tmpDEploidIO.setkStrain(goodProp.size());
  tmpDEploidIO.setInitialPropWasGiven(true);
  tmpDEploidIO.setInitialProp(goodProp);
  tmpDEploidIO.finalProp_ = goodProp;
  tmpDEploidIO.mixture_control_ = mixture_control_;
  tmpDEploidIO.hapParameters() = hapParameters();
  tmpDEploidIO.ibdParameters() = ibdParameters();


  IBDpath tmpIBDpath;
  tmpIBDpath.init(tmpDEploidIO);
  tmpIBDpath.buildPathProbabilityForPainting(goodProp);
  ibdLLK_ = tmpIBDpath.bestPath(goodProp);
  ibdProbsHeader_ = tmpIBDpath.getIBDprobsHeader();
  getIBDprobsIntegrated(tmpIBDpath.getFwdBwd());
  writeIBDpostProb(tmpIBDpath.getFwdBwd(), ibdProbsHeader_);

}

