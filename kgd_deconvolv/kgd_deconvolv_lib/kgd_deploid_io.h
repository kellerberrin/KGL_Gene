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
#include "kgd_ctl_parameter.h"
#include "kgd_random.h"


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
  DEploidIO(const MixtureDataObj& mixture_data);
  ~DEploidIO() = default;

  // Panel.
  std::shared_ptr<Panel> getPanel() { return panel_; }

  // Data
  const std::vector<size_t>& indexOfChromStarts() const { return getMixtureData().indexOfChromStarts(); }
  const std::vector<size_t>& getIndexPosition(size_t index) const { return getMixtureData().getPosition()[index]; }
  const std::vector<std::vector<size_t>>& getPosition() const { return getMixtureData().getPosition(); }
  const std::vector<double>& getPlaf() const { return getMixtureData().getPlaf(); }
  const std::vector<double>& getRefCount() const { return getMixtureData().getRefCount(); }
  const std::vector<double>& getAltCount() const { return getMixtureData().getAltCount(); }
  size_t nLoci() const { return getMixtureData().nLoci(); }
  const MixtureDataObj& getMixtureData() const { return mixture_data_; }
  void setMixtureData(const MixtureDataObj& setTo) { mixture_data_ = setTo; }


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
  const MixtureControlObj& getMixtureControl() const { return mixture_control_; }
  MixtureControlObj& getMixtureControl() { return mixture_control_; }

  // Parameters for the haplotype MCMC
  const HapParameterObj& hapParameters() const { return hap_parameters_; }
  HapParameterObj& hapParameters() { return hap_parameters_; }
  // Parameters for the IBD MCMC
  const IBDParameterObj& ibdParameters() const { return ibd_parameters_; }
  IBDParameterObj& ibdParameters() { return ibd_parameters_; }

  const std::vector<std::vector<double> >& getInitialHap() const { return initialHap_; }
  const std::vector<double>& getInitialProp() const { return  initialProp_; }

  size_t kStrain() const { return k_strain_; }
  void setkStrain(size_t setTo) { k_strain_ = setTo; }


  // Parameter Modifier.
  void writeLastSingleFwdProb(const std::vector<std::vector<double> > &probabilities, size_t chromIndex, size_t strainIndex, bool useIBD);
  void writeMcmcRelated(std::shared_ptr<McmcSample> mcmcSample, bool useIBD = false);
  void setacceptRatio(const double setTo) { acceptRatio_ = setTo; }
  void setmeanThetallks(const double setTo) { meanThetallks_ = setTo; }
  void setmaxLLKs(const double setTo) { maxLLKs_ = setTo; }
  void setmeanllks(const double setTo) { meanllks_ = setTo; }
  void setstdvllks(const double setTo) { stdvllks_ = setTo; }
  void setdicByVar(const double setTo) { dicByVar_ = setTo; }
  void setdicByTheta(const double setTo) { dicByTheta_ = setTo; }

  void setFinalProp(const std::vector<double>& finalProp) { finalProp_ = finalProp; }
  void setInitialProp(const std::vector<double>& initial_prop) { initialProp_ = initial_prop; }
  void setInitialHap(const std::vector<std::vector<double> >& initial_hap) { initialHap_ = initial_hap; }

  // Painting related
  void chromPainting();
  void computeLLKfromInitialHap();
  void paintIBD();
  void getIBDprobsIntegrated(const std::vector<std::vector<double> > &prob);

  // Log
  void wrapUp();

private:

  // Objects.
  std::shared_ptr<Panel> panel_;
  MixtureDataObj mixture_data_;
  MixtureControlObj mixture_control_;
  HapParameterObj hap_parameters_; // Parameters for the haplotype MCMC
  IBDParameterObj ibd_parameters_;   // Parameters for the IBD MCMC

  // Max Number of strains
  size_t k_strain_;

  // Read in input
  std::string plafFileName_;
  std::string refFileName_;
  std::string altFileName_;
  std::string vcfFileName_;
  std::string excludeFileName_;
  std::string initialHapFileName_;
  std::string prefix_;

  std::string panelFileName_;


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
  void readFiles();
  void checkInput();
  void finalize();
  void removeFilesWithSameName();
  std::vector<double> computeExpectedWsafFromInitialHap();
  void computeEffectiveKstrain(std::vector<double> proportion);
  void computeInferredKstrain(std::vector<double> proportion);
  void computeAdjustedEffectiveKstrain();
  void readPanel();

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
