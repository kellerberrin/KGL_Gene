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


#ifndef KGD_PARAM_H
#define KGD_PARAM_H



#include <memory>
#include <vector>
#include <fstream>
#include <stdlib.h>             // strtol, strtod
#include <stdexcept>            // std::invalid_argument
#include <iostream>             // std::cout
#include <sstream>              // std::stringstream
#include "kgd_global.h"
#include "kgd_exceptions.h"
#include "kgd_panel.h"
#include "kgd_vcfReader.h"
#include "kgd_random_generator.h"



namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



class McmcSample;

class UpdateSingleHap;

class UpdatePairHap;

class DEploidIO {
#ifdef UNITTEST
  friend class TestIO;
  friend class TestMcmcMachinery;
  friend class TestIBDpath;
#endif

  friend class McmcMachinery;

  friend class RMcmcSample;

  friend class IBDpath;

public:

  DEploidIO();
  ~DEploidIO() = default;

  // Painting related
  void chromPainting(std::shared_ptr<RandomGenerator> random_generator);

  bool doLsPainting() const { return this->doLsPainting_; }

  bool doIbdPainting() const { return this->doIbdPainting_; }

  bool doComputeLLK() const { return this->doComputeLLK_; }

  void computeLLKfromInitialHap();

  bool useIBD() const { return this->useIBD_; }

  void paintIBD(std::shared_ptr<RandomGenerator> random_generator);

  double ibdLLK_;

  void getIBDprobsIntegrated(std::vector<std::vector<double> > &prob);

  // Log
  void wrapUp();

  size_t randomSeed() const { return randomSeed_; }

private:

  // Object pointers.
  std::shared_ptr<VcfReader> vcfReaderPtr_;
  std::shared_ptr<Panel> panel;
  std::shared_ptr<ExcludeMarker> excludedMarkers;

  double llkFromInitialHap_;

  // Read in input
  std::string plafFileName_;
  std::string refFileName_;
  std::string altFileName_;
  std::string vcfFileName_;
  std::string excludeFileName_;
  std::string initialHapFileName_;
  std::string prefix_;
  size_t randomSeed_;
  bool randomSeedWasGiven_;

  void setrandomSeedWasGiven(const bool random) { this->randomSeedWasGiven_ = random; }


  bool initialPropWasGiven_;
  bool pleaseCheckInitialP_;
  bool initialHapWasGiven_;
  bool kStrainWasManuallySet_;
  bool kStrainWasSetByHap_;
  bool kStrainWasSetByProp_;
  bool useConstRecomb_;
  bool forbidCopyFromSame_;
  size_t kStrain_;
  size_t precision_;
  size_t nMcmcSample_;
  size_t mcmcMachineryRate_;
  double mcmcBurn_;

  bool doUpdateProp_;
  bool doUpdatePair_;
  bool doUpdateSingle_;
  bool doExportPostProb_;
  bool doExportSwitchMissCopy_;
  bool doAllowInbreeding_;
  bool doLsPainting_;
  bool doIbdPainting_;
  bool useIBD_;

  std::vector<double> initialProp_;
  std::vector<double> finalProp_;
  std::vector<std::vector<double> > initialHap_;
  std::vector<std::string> chrom_;
  std::vector<size_t> indexOfChromStarts_;
  std::vector<std::vector<int> > position_;
  std::vector<double> plaf_;
  std::vector<double> refCount_;
  std::vector<double> altCount_;
  size_t nLoci_;


  // Panel related
  bool usePanel_;

  void setUsePanel(const bool setTo) { this->usePanel_ = setTo; }


  bool useVcf_;

  void setUseVcf(const bool useVcf) { this->useVcf_ = useVcf; }

  bool useVcf() const { return this->useVcf_; }

  bool doExportVcf_;

  void setDoExportVcf(const bool exportVcf) { this->doExportVcf_ = exportVcf; }

  bool doExportVcf() const { return this->doExportVcf_; }

  bool compressVcf_;

  void setCompressVcf(const bool compress) { this->compressVcf_ = compress; }

  bool compressVcf() const { return this->compressVcf_; }

  bool doExportRecombProb_;

  void setDoExportRecombProb(const bool exportRecombProb) { this->doExportRecombProb_ = exportRecombProb; }

  bool doExportRecombProb() const { return this->doExportRecombProb_; }

  bool doComputeLLK_;

  void setDoComputeLLK(const bool setTo) { this->doComputeLLK_ = setTo; }

  // Parameters
  double missCopyProb_;
  double averageCentimorganDistance_;// = 15000.0,
  //double Ne_;// = 10.0
  double constRecombProb_;
  double scalingFactor_; // 100.0

  // Diagnostics
  double maxLLKs_;

  void setmaxLLKs(const double setTo) { this->maxLLKs_ = setTo; }

  double meanThetallks_;

  void setmeanThetallks(const double setTo) { this->meanThetallks_ = setTo; }

  double meanllks_;

  void setmeanllks(const double setTo) { this->meanllks_ = setTo; }

  double stdvllks_;

  void setstdvllks(const double setTo) { this->stdvllks_ = setTo; }

  double dicByTheta_;

  void setdicByTheta(const double setTo) { this->dicByTheta_ = setTo; }

  double dicByVar_;

  void setdicByVar(const double setTo) { this->dicByVar_ = setTo; }

  double acceptRatio_;

  void setacceptRatio(const double setTo) { this->acceptRatio_ = setTo; }

  // Output stream
  std::string dEploidGitVersion_;
  std::string compileTime_;
  std::string strExportLLK_;
  std::string strExportHap_;
  std::string strExportVcf_;
  std::string strExportProp_;
  std::string strExportLog_;
  std::string strExportRecombProb_;

  std::string strIbdExportProp_;
  std::string strIbdExportLLK_;
  std::string strIbdExportHap_;
  std::string strIbdExportProbs_;

  std::string strExportSingleFwdProbPrefix_;
  std::string strExportPairFwdProb_;
  std::string strIbdExportSingleFwdProbPrefix_;
  std::string strIbdExportPairFwdProb_;

  std::string strExportExtra_;

  std::ofstream ofstreamExportTmp_;
  std::ofstream ofstreamExportFwdProb_;

  std::string startingTime_;
  std::string endTime_;

  void getTime(bool isStartingTime);


  // Methods
  void init();

  void parse();

  void checkInput();

  void finalize();

  void set_seed(const size_t seed) { this->randomSeed_ = seed; }

  void removeFilesWithSameName();

  std::vector<double> computeExpectedWsafFromInitialHap();


  // Getters and Setters
  void setDoUpdateProp(const bool setTo) { this->doUpdateProp_ = setTo; }

  bool doUpdateProp() const { return this->doUpdateProp_; }

  void setDoUpdateSingle(const bool setTo) { this->doUpdateSingle_ = setTo; }

  bool doUpdateSingle() const { return this->doUpdateSingle_; }

  void setDoUpdatePair(const bool setTo) { this->doUpdatePair_ = setTo; }

  bool doUpdatePair() const { return this->doUpdatePair_; }

  void setDoExportPostProb(const bool setTo) { this->doExportPostProb_ = setTo; }

  bool doExportPostProb() const { return this->doExportPostProb_; }

  void setDoExportSwitchMissCopy(const bool setTo) { this->doExportSwitchMissCopy_ = setTo; }

  bool doExportSwitchMissCopy() const { return this->doExportSwitchMissCopy_; }

  void setDoAllowInbreeding(const bool setTo) { this->doAllowInbreeding_ = setTo; }

  bool doAllowInbreeding() const { return this->doAllowInbreeding_; }

  void setDoLsPainting(const bool setTo) { this->doLsPainting_ = setTo; }

  void setDoIbdPainting(const bool setTo) { this->doIbdPainting_ = setTo; }

  void setUseIBD(const bool setTo) { this->useIBD_ = setTo; }

  bool initialPropWasGiven() const { return initialPropWasGiven_; }

  void setInitialPropWasGiven(const bool setTo) { this->initialPropWasGiven_ = setTo; }

  bool pleaseCheckInitialP() const { return pleaseCheckInitialP_; }

  void setPleaseCheckInitialP(const bool setTo) { this->pleaseCheckInitialP_ = setTo; }

  bool initialHapWasGiven() const { return initialHapWasGiven_; }

  void setInitialHapWasGiven(const bool setTo) { this->initialHapWasGiven_ = setTo; }

  bool randomSeedWasGiven() const { return this->randomSeedWasGiven_; }

  // log and export resutls
  void writeRecombProb(std::shared_ptr<Panel> panel);

  void writeIBDpostProb(std::vector<std::vector<double> > &reshapedProbs, std::vector<std::string> header);

  std::vector<std::string> ibdProbsHeader_;
  std::vector<double> ibdProbsIntegrated_;

  void writeLLK(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);

  void writeProp(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);

  void writeHap(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);

  void writeVcf(std::shared_ptr<McmcSample> mcmcSample);

  void
  writeLastSingleFwdProb(std::vector<std::vector<double> > &probabilities, size_t chromIndex, size_t strainIndex, bool useIBD);

  void writeLastPairFwdProb(UpdatePairHap &updatePair, size_t chromIndex);

  void writeLog(std::ostream *writeTo);

  void writeEventCount();

  std::vector<double> IBDpathChangeAt_;
  std::vector<double> finalIBDpathChangeAt_;

  std::vector<double> siteOfTwoSwitchOne_;
  std::vector<double> siteOfTwoMissCopyOne_;
  std::vector<double> siteOfTwoSwitchTwo_;
  std::vector<double> siteOfTwoMissCopyTwo_;
  std::vector<double> siteOfOneSwitchOne_;
  std::vector<double> siteOfOneMissCopyOne_;

  std::vector<double> finalSiteOfTwoSwitchOne_;
  std::vector<double> finalSiteOfTwoMissCopyOne_;
  std::vector<double> finalSiteOfTwoSwitchTwo_;
  std::vector<double> finalSiteOfTwoMissCopyTwo_;
  std::vector<double> finalSiteOfOneSwitchOne_;
  std::vector<double> finalSiteOfOneMissCopyOne_;


  void writeMcmcRelated(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);

  void readPanel();

  // Panel related
  bool usePanel() const { return usePanel_; }

  std::string panelFileName_;
  double parameterG_;

  void setParameterG(const double setTo) { this->parameterG_ = setTo; }

  double parameterG() const { return this->parameterG_; }

  double parameterSigma_;

  void setParameterSigma(const double setTo) { this->parameterSigma_ = setTo; }

  double parameterSigma() const { return this->parameterSigma_; }

  double ibdSigma_;

  void setIBDSigma(const double setTo) { this->ibdSigma_ = setTo; }

  double ibdSigma() const { return this->ibdSigma_; }

  void setNLoci(const size_t setTo) { this->nLoci_ = setTo; }

  size_t nLoci() const { return this->nLoci_; }

  void setKstrain(const size_t setTo) { this->kStrain_ = setTo; }

  size_t kStrain() const { return this->kStrain_; }

  void setKStrainWasManuallySet(const size_t setTo) { this->kStrainWasManuallySet_ = setTo; }

  bool kStrainWasSetByHap() const { return this->kStrainWasSetByHap_; }

  void setKStrainWasSetByHap(const size_t setTo) { this->kStrainWasSetByHap_ = setTo; }

  bool kStrainWasManuallySet() const { return this->kStrainWasManuallySet_; }

  void setKStrainWasSetByProp(const size_t setTo) { this->kStrainWasSetByProp_ = setTo; }

  bool kStrainWasSetByProp() const { return this->kStrainWasSetByProp_; }

  size_t nMcmcSample() const { return this->nMcmcSample_; }

  double averageCentimorganDistance() const { return this->averageCentimorganDistance_; }

  //double Ne() const { return this->Ne_; }
  double scalingFactor() const { return this->scalingFactor_; }

  void setScalingFactor(const double setTo) { this->scalingFactor_ = setTo; }

  double constRecombProb() const { return this->constRecombProb_; }

  bool useConstRecomb() const { return this->useConstRecomb_; }

  bool excludeSites_;

  bool excludeSites() const { return this->excludeSites_; }

  void setExcludeSites(const size_t exclude) { this->excludeSites_ = exclude; }

  bool forbidCopyFromSame() const { return this->forbidCopyFromSame_; }

  void setForbidCopyFromSame(const bool forbid) { this->forbidCopyFromSame_ = forbid; }

  double effectiveKstrain_;

  void computeEffectiveKstrain(std::vector<double> proportion);

  int inferredKstrain_;

  void computeInferredKstrain(std::vector<double> proportion);

  double adjustedEffectiveKstrain_;

  void computeAdjustedEffectiveKstrain();

};



}   // organization level namespace
}   // project level namespace



#endif
