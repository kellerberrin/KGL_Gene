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


#ifndef KGD_PARAM_H
#define KGD_PARAM_H



#include <memory>
#include <vector>
#include <fstream>
#include <cstdlib>             // strtol, strtod
#include <stdexcept>            // std::invalid_argument
#include <iostream>             // std::cout
#include <sstream>              // std::stringstream
#include "kgd_global.h"
#include "kgd_exceptions.h"
#include "kgd_panel.h"
#include "kgd_vcf_reader.h"
#include "kgd_ctl_data.h"
#include "kgd_ctl_function.h"
#include "kgd_random_generator.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


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

  // Panel.
  std::shared_ptr<Panel> getPanel() { return panel_; }

  // Data
  const std::vector<size_t>& indexOfChromStarts() const { return getMixtureData().indexOfChromStarts(); }
  const std::vector<int>& getIndexPosition(size_t index) const { return getMixtureData().getPosition()[index]; }
  const std::vector<std::vector<int>>& getPosition() const { return getMixtureData().getPosition(); }
  const std::vector<double>& getPlaf() const { return getMixtureData().getPlaf(); }
  const std::vector<double>& getRefCount() const { return getMixtureData().getRefCount(); }
  const std::vector<double>& getAltCount() const { return getMixtureData().getAltCount(); }
  size_t nLoci() const { return getMixtureData().nLoci(); }
  const MixtureDataObj& getMixtureData() const { return mixture_data; }

  // Set Control
  void setInitialPropWasGiven(const bool setTo) { getMixtureControl().setInitialPropWasGiven(setTo); }
  void setDoUpdateProp(const bool setTo) { getMixtureControl().setDoUpdateProp(setTo); }
  void setInitialHapWasGiven(const bool setTo) { getMixtureControl().setInitialHapWasGiven(setTo); }
  // Get Control
  bool doExportPostProb() const { return getMixtureControl().doExportPostProb(); }
  bool doAllowInbreeding() const { return getMixtureControl().doAllowInbreeding(); }
  bool initialPropWasGiven() const { return getMixtureControl().initialPropWasGiven(); }
  bool initialHapWasGiven() const { return getMixtureControl().initialHapWasGiven(); }
  bool doUpdateProp() const { return getMixtureControl().doUpdateProp(); }
  bool doUpdateSingle() const { return getMixtureControl().doUpdateSingle(); }
  bool doUpdatePair() const { return getMixtureControl().doUpdatePair(); }
  bool forbidCopyFromSame() const { return getMixtureControl().forbidCopyFromSame(); }
  bool useConstRecomb() const { return getMixtureControl().useConstRecomb(); }
  bool usePanel() const { return getMixtureControl().usePanel(); }
  bool useIBD() const { return getMixtureControl().useIBD(); }
  const MixtureControlObj& getMixtureControl() const { return mixture_control; }
  MixtureControlObj& getMixtureControl() { return mixture_control; }

  // Get results.

  // Get Parameters
  size_t getMcmcSample() const { return nMcmcSample_; }
  size_t kStrain() const { return kStrain_; }
  size_t getMcmcMachineryRate() const { return mcmcMachineryRate_; }
  double getMcmcBurn() const { return mcmcBurn_; }
  double scalingFactor() const { return scalingFactor_; }
  double getMissCopyProb() const { return missCopyProb_; }
  double ibdSigma() const { return ibdSigma_; }
  double parameterSigma() const { return parameterSigma_; }
  const std::vector<std::vector<double> >& getInitialHap() const { return initialHap_; }
  const std::vector<double>& getInitialProp() const { return  initialProp_; }
  double averageCentimorganDistance() const { return averageCentimorganDistance_; }
  double parameterG() const { return parameterG_; }
  double constRecombProb() const { return constRecombProb_; }
  size_t randomSeed() const { return randomSeed_; }

  // Parameter Modifier.
  void writeLastSingleFwdProb(const std::vector<std::vector<double> > &probabilities, size_t chromIndex, size_t strainIndex, bool useIBD);
  void setFinalProp(const std::vector<double>& finalProp) { finalProp_ = finalProp; }
  void writeMcmcRelated(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);
  void setInitialProp(const std::vector<double>& initial_prop) { initialProp_ = initial_prop; }
  void setInitialHap(const std::vector<std::vector<double> >& initial_hap) { initialHap_ = initial_hap; }
  void setacceptRatio(const double setTo) { acceptRatio_ = setTo; }
  void setmeanThetallks(const double setTo) { meanThetallks_ = setTo; }
  void setmaxLLKs(const double setTo) { maxLLKs_ = setTo; }
  void setmeanllks(const double setTo) { meanllks_ = setTo; }
  void setstdvllks(const double setTo) { stdvllks_ = setTo; }
  void setdicByVar(const double setTo) { dicByVar_ = setTo; }
  void setdicByTheta(const double setTo) { dicByTheta_ = setTo; }


  // Painting related
  void chromPainting(std::shared_ptr<RandomGenerator> random_generator);
  void computeLLKfromInitialHap();
  void paintIBD(std::shared_ptr<RandomGenerator> random_generator);
  void getIBDprobsIntegrated(const std::vector<std::vector<double> > &prob);

  // Log
  void wrapUp();

private:

  // Objects.
  std::shared_ptr<Panel> panel_;
  MixtureDataObj mixture_data;
  MixtureControlObj mixture_control;

  // Read in input
  std::string plafFileName_;
  std::string refFileName_;
  std::string altFileName_;
  std::string vcfFileName_;
  std::string excludeFileName_;
  std::string initialHapFileName_;
  std::string prefix_;

  std::string panelFileName_;


  size_t nMcmcSample_;    // Number of MCMC samples (default value 800).
  size_t mcmcMachineryRate_;  // MCMC sample rate (default value 5).
  double mcmcBurn_;  // MCMC burn rate (default value 0.5).

  size_t kStrain_;
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
  double ibdLLK_;
  double llkFromInitialHap_;

  std::vector<std::string> ibdProbsHeader_;
  std::vector<double> ibdProbsIntegrated_;

  std::vector<double> initialProp_;
  std::vector<double> finalProp_;
  std::vector<std::vector<double> > initialHap_;

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
  void removeFilesWithSameName();
  std::vector<double> computeExpectedWsafFromInitialHap();
  void computeEffectiveKstrain(std::vector<double> proportion);
  void computeInferredKstrain(std::vector<double> proportion);
  void computeAdjustedEffectiveKstrain();
  void readPanel();

  void set_seed(const size_t seed) { randomSeed_ = seed; }
  void setParameterG(const double setTo) { parameterG_ = setTo; }
  void setParameterSigma(const double setTo) { parameterSigma_ = setTo; }
  void setIBDSigma(const double setTo) { ibdSigma_ = setTo; }
  void setKstrain(const size_t setTo) { kStrain_ = setTo; }
  void setScalingFactor(const double setTo) { scalingFactor_ = setTo; }

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
