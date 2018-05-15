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
#include "kgd_vcf_reader.h"
#include "kgd_random_generator.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace


// Forward declarations.
class McmcSample;
class UpdateSingleHap;
class UpdatePairHap;

class DEploidIO {
#ifdef UNITTEST
  friend class TestIO;
  friend class TestMcmcMachinery;
  friend class TestIBDpath;
#endif


public:

  DEploidIO();
  ~DEploidIO() = default;



  // Access Routines.
  std::shared_ptr<Panel> getPanel() { return panel_; }

  size_t randomSeed() const { return randomSeed_; }
  bool doExportPostProb() const { return doExportPostProb_; }
  bool doAllowInbreeding() const { return doAllowInbreeding_; }
  bool initialPropWasGiven() const { return initialPropWasGiven_; }
  bool initialHapWasGiven() const { return initialHapWasGiven_; }
  bool doUpdateProp() const { return doUpdateProp_; }
  bool doUpdateSingle() const { return doUpdateSingle_; }
  bool doUpdatePair() const { return doUpdatePair_; }
  bool forbidCopyFromSame() const { return forbidCopyFromSame_; }
  bool useConstRecomb() const { return useConstRecomb_; }

  const std::vector<size_t>& indexOfChromStarts() const { return indexOfChromStarts_; }
  const std::vector<int>& getIndexPosition(size_t index) const { return position_[index]; }
  const std::vector<std::vector<int>>& getPosition() const { return position_; }
  const std::vector<double>& getPlaf() const { return plaf_; }
  const std::vector<double>& getRefCount() const { return refCount_; }
  const std::vector<double>& getAltCount() const { return altCount_; }

  double scalingFactor() const { return scalingFactor_; }
  double getMissCopyProb() const { return missCopyProb_; }
  size_t getMcmcSample() const { return nMcmcSample_; }
  size_t getMcmcMachineryRate() const { return mcmcMachineryRate_; }
  double getMcmcBurn() const { return mcmcBurn_; }
  size_t kStrain() const { return kStrain_; }
  double ibdSigma() const { return ibdSigma_; }
  double parameterSigma() const { return parameterSigma_; }
  const std::vector<std::vector<double> >& getInitialHap() const { return initialHap_; }
  const std::vector<double>& getInitialProp() const { return  initialProp_; }
  size_t nLoci() const { return nLoci_; }
  double averageCentimorganDistance() const { return averageCentimorganDistance_; }
  double parameterG() const { return parameterG_; }
  double constRecombProb() const { return constRecombProb_; }

  // Modifiers.
  void writeLastSingleFwdProb(const std::vector<std::vector<double> > &probabilities, size_t chromIndex, size_t strainIndex, bool useIBD);
  void setFinalProp(const std::vector<double>& finalProp) { finalProp_ = finalProp; }
  void writeMcmcRelated(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);
  void setInitialProp(const std::vector<double>& initial_prop) { initialProp_ = initial_prop; }
  void setInitialPropWasGiven(const bool setTo) { initialPropWasGiven_ = setTo; }
  void setDoUpdateProp(const bool setTo) { doUpdateProp_ = setTo; }
  void setInitialHap(const std::vector<std::vector<double> >& initial_hap) { initialHap_ = initial_hap; }
  void setInitialHapWasGiven(const bool setTo) { initialHapWasGiven_ = setTo; }
  void setacceptRatio(const double setTo) { acceptRatio_ = setTo; }
  void setmeanThetallks(const double setTo) { meanThetallks_ = setTo; }
  void setmaxLLKs(const double setTo) { maxLLKs_ = setTo; }
  void setmeanllks(const double setTo) { meanllks_ = setTo; }
  void setstdvllks(const double setTo) { stdvllks_ = setTo; }
  void setdicByVar(const double setTo) { dicByVar_ = setTo; }
  void setdicByTheta(const double setTo) { dicByTheta_ = setTo; }

  // Painting related
  void chromPainting(std::shared_ptr<RandomGenerator> random_generator);
  bool doLsPainting() const { return doLsPainting_; }
  bool doIbdPainting() const { return doIbdPainting_; }
  bool doComputeLLK() const { return doComputeLLK_; }
  void computeLLKfromInitialHap();
  bool useIBD() const { return useIBD_; }
  void paintIBD(std::shared_ptr<RandomGenerator> random_generator);
  void getIBDprobsIntegrated(const std::vector<std::vector<double> > &prob);

  // Log
  void wrapUp();

private:

  // Object pointers.
  std::shared_ptr<VcfReader> vcfReaderPtr_;
  std::shared_ptr<Panel> panel_;
  std::shared_ptr<ExcludeMarker> excludedMarkers_;

  // Read in input
  std::string plafFileName_;
  std::string refFileName_;
  std::string altFileName_;
  std::string vcfFileName_;
  std::string excludeFileName_;
  std::string initialHapFileName_;
  std::string prefix_;

  std::vector<std::string> ibdProbsHeader_;
  std::vector<double> ibdProbsIntegrated_;
  std::string panelFileName_;

  bool randomSeedWasGiven_;
  bool initialPropWasGiven_;
  bool pleaseCheckInitialP_;
  bool initialHapWasGiven_;
  bool kStrainWasManuallySet_;
  bool kStrainWasSetByHap_;
  bool kStrainWasSetByProp_;
  bool useConstRecomb_;
  bool forbidCopyFromSame_;
  bool doUpdateProp_;
  bool doUpdatePair_;
  bool doUpdateSingle_;
  bool doExportPostProb_;
  bool doExportSwitchMissCopy_;
  bool doAllowInbreeding_;
  bool doLsPainting_;
  bool doIbdPainting_;
  bool useIBD_;
  bool usePanel_;  // Panel related
  bool useVcf_;
  bool doExportVcf_;
  bool compressVcf_;
  bool doExportRecombProb_;
  bool doComputeLLK_;
  bool excludeSites_;


  size_t kStrain_;
  size_t precision_;
  size_t nMcmcSample_;
  size_t mcmcMachineryRate_;
  double mcmcBurn_;
  size_t nLoci_;
  double ibdLLK_;
  double llkFromInitialHap_;
  size_t randomSeed_;

  // Parameters
  double missCopyProb_;
  double averageCentimorganDistance_;// = 15000.0,
  double constRecombProb_;
  double scalingFactor_; // 100.0

  double parameterG_;
  double parameterSigma_;
  double ibdSigma_;

  // Diagnostics
  double maxLLKs_;
  double meanThetallks_;
  double meanllks_;
  double stdvllks_;
  double dicByTheta_;
  double dicByVar_;
  double acceptRatio_;

  // Computed values.
  double effectiveKstrain_;
  int inferredKstrain_;
  double adjustedEffectiveKstrain_;

  std::vector<double> initialProp_;
  std::vector<double> finalProp_;
  std::vector<std::vector<double> > initialHap_;
  std::vector<std::string> chrom_;
  std::vector<size_t> indexOfChromStarts_;
  std::vector<std::vector<int> > position_;
  std::vector<double> plaf_;
  std::vector<double> refCount_;
  std::vector<double> altCount_;


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



  // Methods
  void getTime(bool isStartingTime);
  void init();
  void parse();
  void checkInput();
  void finalize();
  void set_seed(const size_t seed) { randomSeed_ = seed; }
  void removeFilesWithSameName();
  std::vector<double> computeExpectedWsafFromInitialHap();
  void computeEffectiveKstrain(std::vector<double> proportion);
  void computeInferredKstrain(std::vector<double> proportion);
  void computeAdjustedEffectiveKstrain();
  void readPanel();

  // Getters and Setters
  void setDoUpdateSingle(const bool setTo) { doUpdateSingle_ = setTo; }
  void setDoUpdatePair(const bool setTo) { doUpdatePair_ = setTo; }
  void setDoExportPostProb(const bool setTo) { doExportPostProb_ = setTo; }
  void setDoExportSwitchMissCopy(const bool setTo) { doExportSwitchMissCopy_ = setTo; }
  bool doExportSwitchMissCopy() const { return doExportSwitchMissCopy_; }
  void setDoAllowInbreeding(const bool setTo) { doAllowInbreeding_ = setTo; }
  void setDoLsPainting(const bool setTo) { doLsPainting_ = setTo; }
  void setDoIbdPainting(const bool setTo) { doIbdPainting_ = setTo; }
  void setUseIBD(const bool setTo) { useIBD_ = setTo; }
  bool pleaseCheckInitialP() const { return pleaseCheckInitialP_; }
  void setPleaseCheckInitialP(const bool setTo) { pleaseCheckInitialP_ = setTo; }
  bool randomSeedWasGiven() const { return randomSeedWasGiven_; }
  void setUsePanel(const bool setTo) { usePanel_ = setTo; }
  void setUseVcf(const bool useVcf) { useVcf_ = useVcf; }
  bool useVcf() const { return useVcf_; }
  void setDoExportVcf(const bool exportVcf) { doExportVcf_ = exportVcf; }
  bool doExportVcf() const { return doExportVcf_; }
  void setCompressVcf(const bool compress) { compressVcf_ = compress; }
  bool compressVcf() const { return compressVcf_; }
  void setDoExportRecombProb(const bool exportRecombProb) { doExportRecombProb_ = exportRecombProb; }
  bool doExportRecombProb() const { return doExportRecombProb_; }
  void setDoComputeLLK(const bool setTo) { doComputeLLK_ = setTo; }
  void setrandomSeedWasGiven(const bool random) { randomSeedWasGiven_ = random; }
  void setParameterG(const double setTo) { parameterG_ = setTo; }
  void setParameterSigma(const double setTo) { parameterSigma_ = setTo; }
  void setIBDSigma(const double setTo) { ibdSigma_ = setTo; }
  void setNLoci(const size_t setTo) { nLoci_ = setTo; }
  void setKstrain(const size_t setTo) { kStrain_ = setTo; }
  void setKStrainWasManuallySet(const size_t setTo) { kStrainWasManuallySet_ = setTo; }
  bool kStrainWasSetByHap() const { return kStrainWasSetByHap_; }
  void setKStrainWasSetByHap(const size_t setTo) { kStrainWasSetByHap_ = setTo; }
  bool kStrainWasManuallySet() const { return kStrainWasManuallySet_; }
  void setKStrainWasSetByProp(const size_t setTo) { kStrainWasSetByProp_ = setTo; }
  bool kStrainWasSetByProp() const { return kStrainWasSetByProp_; }
  void setScalingFactor(const double setTo) { scalingFactor_ = setTo; }
  bool excludeSites() const { return excludeSites_; }
  void setExcludeSites(const size_t exclude) { excludeSites_ = exclude; }
  void setForbidCopyFromSame(const bool forbid) { forbidCopyFromSame_ = forbid; }
  bool usePanel() const { return usePanel_; }

  // log and export results
  void writeRecombProb(std::shared_ptr<Panel> panel);
  void writeIBDpostProb(const std::vector<std::vector<double> > &reshapedProbs, std::vector<std::string> header);
  void writeLLK(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);
  void writeProp(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);
  void writeHap(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);
  void writeVcf(std::shared_ptr<McmcSample> mcmcSample);
  void writeLastPairFwdProb(UpdatePairHap &updatePair, size_t chromIndex);
  void writeLog(std::ostream *writeTo);
  void writeEventCount();


};



}   // organization level namespace
}   // project level namespace



#endif
