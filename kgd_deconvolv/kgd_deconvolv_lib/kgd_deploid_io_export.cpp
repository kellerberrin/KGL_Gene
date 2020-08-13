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

#include <iomanip>      // std::setw
#include "kgd_deploid_io.h"
#include "kgd_deconvolv_app.h"


namespace kgd = kellerberrin::deconvolv;


void kgd::DEploidIO::wrapUp() {

  writeRecombProb(panel_);

  // Get End time before writing the log
  getTime(false);

  writeLog(&std::cout);

  ofstreamExportTmp_.open(strExportLog_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  writeLog(&ofstreamExportTmp_);

  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeRecombProb(std::shared_ptr<Panel> panel) {

  if (!getMixtureControl().doExportRecombProb()) return;

  if (panel) {

    ofstreamExportTmp_.open(strExportRecombProb_.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    ofstreamExportTmp_ << "p.recomb" << "\t"
                      << "p.each" << "\t"
                      << "p.no.recomb" << "\t"
                      << "p.rec.rec" << "\t"
                      << "p.rec.norec" << "\t"
                      << "p.norec.norec" << "\n";

    for (size_t i = 0; i < panel->getRec().size(); i++) {

      ofstreamExportTmp_ << panel->getRecIndex(i) << "\t"
                        << panel->getRecEachHapIndex(i) << "\t"
                        << panel->getNoRecIndex(i) << "\t"
                        << panel->getRecRecIndex(i) << "\t"
                        << panel->getRecNoRecIndex(i) << "\t"
                        << panel->getNoRecNoRecIndex(i) << "\n";
    }

    ofstreamExportTmp_.close();

  }

}


void kgd::DEploidIO::writeLog(std::ostream *writeTo) {

  size_t nHash = 30 + std::string(Deconvolv::VERSION).size();

  for (size_t i = 0; i < nHash; i++) {

    (*writeTo) << "#";

  }

  (*writeTo) << "\n";
  (*writeTo) << "#        kgd_deconvolv " << std::setw(10) << Deconvolv::VERSION << " log        #\n";

  for (size_t i = 0; i < nHash; i++) {

    (*writeTo) << "#";

  }

  (*writeTo) << "\n";
  (*writeTo) << "Program was compiled on: " << compileTime_ << std::endl;
  (*writeTo) << "kgd_deconvolv version: " << dEploidGitVersion_ << std::endl;
  (*writeTo) << "\n";
  (*writeTo) << "Input data: \n";

  if (panelFileName_.size() > 0) {

    (*writeTo) << std::setw(12) << "Panel: " << panelFileName_ << "\n";

  }

  (*writeTo) << std::setw(12) << "PLAF: " << plafFileName_ << "\n";

  if (getMixtureControl().useVcf()) (*writeTo) << std::setw(12) << "VCF: " << vcfFileName_ << "\n";

  if (refFileName_.size() > 0) (*writeTo) << std::setw(12) << "REF count: " << refFileName_ << "\n";

  if (altFileName_.size() > 0) (*writeTo) << std::setw(12) << "ALT count: " << altFileName_ << "\n";

  if (getMixtureControl().excludeSites()) { (*writeTo) << std::setw(12) << "Exclude: " << excludeFileName_ << "\n"; }

  (*writeTo) << "\n";

  if ((not getMixtureControl().doLsPainting()) && (not getMixtureControl().doIbdPainting())) {

    (*writeTo) << "MCMC parameters: " << "\n";
    (*writeTo) << std::setw(19) << " MCMC burn: " << hapParameters().McmcBurn() << "\n";
    (*writeTo) << std::setw(19) << " MCMC sample: " << hapParameters().McmcSample() << "\n";
    (*writeTo) << std::setw(19) << " MCMC sample rate: " << hapParameters().McmcMachineryRate() << "\n";
    (*writeTo) << std::setw(19) << " Random seed: " << hapParameters().randomSeed() << "\n";

    if (useIBD()) {

      (*writeTo) << std::setw(19) << "  IBD Method used: YES" << "\n";

    }

    (*writeTo) << std::setw(19) << " Update Prop: " << (doUpdateProp() ? "YES" : "NO") << "\n";
    (*writeTo) << std::setw(19) << " Update Single: " << (doUpdateSingle() ? "YES" : "NO") << "\n";
    (*writeTo) << std::setw(19) << " Update Pair: " << (doUpdatePair() ? "YES" : "NO") << "\n";
    (*writeTo) << "\n";

  }

  (*writeTo) << "Other parameters:" << "\n";

  if (getMixtureControl().forbidCopyFromSame()) {

    (*writeTo) << " Update pair haplotypes move forbid copying from the same strain!!! \n";

  }

  (*writeTo) << std::setw(20) << " Miss copy prob: " << hapParameters().getMissCopyProb() << "\n";
  (*writeTo) << std::setw(20) << " Avrg Cent Morgan: " << hapParameters().averageCentimorganDistance() << "\n";
  (*writeTo) << std::setw(20) << " G: " << hapParameters().parameterG() << "\n";

  if (useIBD()) {

    (*writeTo) << std::setw(20) << " IBD sigma: " << ibdParameters().proposalSigma() << "\n";

  } else {

    (*writeTo) << std::setw(20) << " sigma: " << hapParameters().proposalSigma() << "\n";

  }

  (*writeTo) << std::setw(20) << " ScalingFactor: " << hapParameters().proposalUpdateScaling() << "\n";

  if (initialPropWasGiven()) {

    (*writeTo) << std::setw(20) << " Initial prob: ";

    for (size_t i = 0; i < initialProp_.size(); i++) {

      (*writeTo) << initialProp_[i] << ((i != (kStrain() - 1)) ? " " : "\n");

    }

  }

  (*writeTo) << "\n";

  if ((not getMixtureControl().doLsPainting())
      and (not getMixtureControl().doIbdPainting())
      and (not getMixtureControl().doComputeLLK())) {

    (*writeTo) << "MCMC diagnostic:" << "\n";
    (*writeTo) << std::setw(19) << " Accept_ratio: " << acceptRatio_ << "\n";
    (*writeTo) << std::setw(19) << " Max_llks: " << maxLLKs_ << "\n";
    (*writeTo) << std::setw(19) << " Final_theta_llks: " << meanThetallks_ << "\n";
    (*writeTo) << std::setw(19) << " Mean_llks: " << meanllks_ << "\n";
    (*writeTo) << std::setw(19) << " Stdv_llks: " << stdvllks_ << "\n";
    (*writeTo) << std::setw(19) << " DIC_by_Dtheta: " << dicByTheta_ << "\n";
    (*writeTo) << std::setw(19) << " DIC_by_varD: " << dicByVar_ << "\n";
    (*writeTo) << "\n";

  }

  (*writeTo) << "Run time:\n";
  (*writeTo) << std::setw(14) << "Start at: " << startingTime_;
  (*writeTo) << std::setw(14) << "End at: " << endTime_;
  (*writeTo) << "\n";

  if (getMixtureControl().doComputeLLK()) {

    (*writeTo) << "Input likelihood: " << llkFromInitialHap_;
    (*writeTo) << "\n";

  } else {

    (*writeTo) << "Output saved to:\n";

    if (getMixtureControl().doLsPainting()) {

      for (size_t i = 0; i < kStrain(); i++) {

        (*writeTo) << "Posterior probability of strain " << i << ": " << strExportSingleFwdProbPrefix_ << i << std::endl;

      }

    } else if (getMixtureControl().doIbdPainting()) {

      if (ibdProbsIntegrated_.size() > 1) {

        (*writeTo) << std::setw(14) << "IBD probs: " << strIbdExportProbs_ << "\n\n";
        (*writeTo) << " IBD probabilities:\n";

        for (size_t stateI = 0; stateI < ibdProbsHeader_.size(); stateI++) {

          (*writeTo) << std::setw(14) << ibdProbsHeader_[stateI] << ": " << ibdProbsIntegrated_[stateI] << "\n";

        }

      }

    } else {

      (*writeTo) << std::setw(14) << "Likelihood: " << strExportLLK_ << "\n";
      (*writeTo) << std::setw(14) << "Proportions: " << strExportProp_ << "\n";
      (*writeTo) << std::setw(14) << "Haplotypes: " << strExportHap_ << "\n";

      if (getMixtureControl().doExportVcf()) { (*writeTo) << std::setw(14) << "Vcf: " << strExportVcf_ << "\n"; }

      if (useIBD()) {

        (*writeTo) << " IBD method output saved to:\n";
        (*writeTo) << std::setw(14) << "Likelihood: " << strIbdExportLLK_ << "\n";
        (*writeTo) << std::setw(14) << "Proportions: " << strIbdExportProp_ << "\n";
        (*writeTo) << std::setw(14) << "Haplotypes: " << strIbdExportHap_ << "\n";

      }

      if (ibdProbsIntegrated_.size() > 1) {

        (*writeTo) << std::setw(14) << "IBD probs: " << strIbdExportProbs_ << "\n\n";
        (*writeTo) << " IBD probabilities:\n";

        for (size_t stateI = 0; stateI < ibdProbsHeader_.size(); stateI++) {

          (*writeTo) << std::setw(14) << ibdProbsHeader_[stateI] << ": " << ibdProbsIntegrated_[stateI] << "\n";

        }

      }

    }

    (*writeTo) << "\n";
    (*writeTo) << " IBD best path llk: " << ibdLLK_ << "\n\n";

    computeEffectiveKstrain(finalProp_);
    (*writeTo) << "         Effective_K: " << effectiveKstrain_ << "\n";
    computeInferredKstrain(finalProp_);
    (*writeTo) << "          Inferred_K: " << inferredKstrain_ << "\n";
    computeAdjustedEffectiveKstrain();
    (*writeTo) << "Adjusted_effective_K: " << adjustedEffectiveKstrain_ << "\n";

  }

  (*writeTo) << "\n";
  (*writeTo) << "Proportions:\n";

  for (size_t ii = 0; ii < finalProp_.size(); ii++) {

    (*writeTo) << std::setw(10) << finalProp_[ii];
    (*writeTo) << ((ii < (finalProp_.size() - 1)) ? "\t" : "\n");

  }

}


void kgd::DEploidIO::writeEventCount() {

  ofstreamExportTmp_.open(strExportExtra_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  // HEADER
  ofstreamExportTmp_ << "CHROM" << "\t"
                    << "POS" << "\t"
                    << "IBD_path_change_at_" << "\t"
                    << "finalIBDpathChangeAt" << "\t"

                    << "siteOfTwoSwitchOne_" << "\t"
                    << "finalSiteOfTwoSwitchOne" << "\t"

                    << "siteOfTwoMissCopyOne_" << "\t"
                    << "finalSiteOfTwoMissCopyOne" << "\t"

                    << "siteOfTwoSwitchTwo_" << "\t"
                    << "finalSiteOfTwoSwitchTwo" << "\t"

                    << "siteOfTwoMissCopyTwo_" << "\t"
                    << "finalSiteOfTwoMissCopyTwo" << "\t"

                    << "siteOfOneSwitchOne" << "\t"
                    << "finalSiteOfOneSwitchOne" << "\t"

                    << "siteOfOneMissCopyOne" << "\t"
                    << "finalSiteOfOneMissCopyOne" << std::endl;

  size_t siteIndex = 0;

  for (size_t chromI = 0; chromI < getMixtureData().getChrom().size(); chromI++) {

    for (size_t posI = 0; posI < getMixtureData().getPosition()[chromI].size(); posI++) {

      ofstreamExportTmp_ << getMixtureData().getChrom()[chromI] << "\t"
                        << static_cast<int>(getMixtureData().getPosition()[chromI][posI]) << "\t"

                        << IBDpathChangeAt_[siteIndex] << "\t"
                        << finalIBDpathChangeAt_[siteIndex] << "\t"

                        << siteOfTwoSwitchOne_[siteIndex] << "\t"
                        << finalSiteOfTwoSwitchOne_[siteIndex] << "\t"

                        << siteOfTwoMissCopyOne_[siteIndex] << "\t"
                        << finalSiteOfTwoMissCopyOne_[siteIndex] << "\t"

                        << siteOfTwoSwitchTwo_[siteIndex] << "\t"
                        << finalSiteOfTwoSwitchTwo_[siteIndex] << "\t"

                        << siteOfTwoMissCopyTwo_[siteIndex] << "\t"
                        << finalSiteOfTwoMissCopyTwo_[siteIndex] << "\t"

                        << siteOfOneSwitchOne_[siteIndex] << "\t"
                        << finalSiteOfOneSwitchOne_[siteIndex] << "\t"

                        << siteOfOneMissCopyOne_[siteIndex] << "\t"
                        << finalSiteOfOneMissCopyOne_[siteIndex] << std::endl;

      siteIndex++;

    }

  }

  assert(siteIndex == IBDpathChangeAt_.size());

  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeIBDpostProb(const std::vector<std::vector<double> > &reshapedProbs, std::vector<std::string> header) {

  std::ostream *writeTo;

#ifdef UNITTEST
  writeTo = &std::cout;
#endif

#ifndef UNITTEST
  ofstreamExportTmp_.open(strIbdExportProbs_.c_str(), std::ios::out | std::ios::app | std::ios::binary);
  writeTo = &ofstreamExportTmp_;
#endif

  (*writeTo) << "CHROM" << "\t" << "POS" << "\t";

  for (std::string tmp : header) {

    (*writeTo) << tmp << ((tmp != header[header.size() - 1]) ? "\t" : "\n");

  }

  size_t siteIndex = 0;

  for (size_t chromIndex = 0; chromIndex < getMixtureData().getPosition().size(); chromIndex++) {

    for (size_t posI = 0; posI < getMixtureData().getPosition()[chromIndex].size(); posI++) {

      (*writeTo) << getMixtureData().getChrom()[chromIndex] << "\t" << static_cast<int>(getMixtureData().getPosition()[chromIndex][posI]) << "\t";

      for (size_t ij = 0; ij < reshapedProbs[siteIndex].size(); ij++) {

        (*writeTo) << reshapedProbs[siteIndex][ij] << "\t";

      }

      (*writeTo) << std::endl;
      siteIndex++;

    }

  }
  assert(siteIndex == nLoci());

#ifndef UNITTEST
  ofstreamExportTmp_.close();
#endif
}


void kgd::DEploidIO::writeLastSingleFwdProb(const std::vector<std::vector<double> > &probabilities,
                                            size_t chromIndex,
                                            size_t strainIndex,
                                            bool useIBD) {

  if (probabilities.size() == 0) {

    return;

  }

  size_t panelSize = probabilities[0].size();

  std::string strExportFwdProb = ((useIBD == true) ? strIbdExportSingleFwdProbPrefix_ : strExportSingleFwdProbPrefix_) + std::to_string(strainIndex);

  ofstreamExportFwdProb_.open(strExportFwdProb.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  if (chromIndex == 0) { // Print header

    ofstreamExportFwdProb_ << "CHROM" << "\t" << "POS" << "\t";

    for (size_t ii = 0; ii < probabilities[0].size(); ii++) {

      if (doAllowInbreeding() == true) {

        if (ii <= (panelSize - kStrain())) {

          ofstreamExportFwdProb_ << "P" << (ii + 1);

        } else {

          ofstreamExportFwdProb_ << "I" << (ii) - (panelSize - kStrain());

        }

      } else {

        ofstreamExportFwdProb_ << (ii + 1);

      }

      ofstreamExportFwdProb_ << ((ii < (panelSize - 1)) ? "\t" : "\n");

    }

  }

  size_t siteIndex = 0;

  for (size_t posI = 0; posI < getMixtureData().getPosition()[chromIndex].size(); posI++) {

    ofstreamExportFwdProb_ << getMixtureData().getChrom()[chromIndex] << "\t" << static_cast<int>(getMixtureData().getPosition()[chromIndex][posI]) << "\t";

    for (size_t ii = 0; ii < probabilities[siteIndex].size(); ii++) {

      ofstreamExportFwdProb_ << probabilities[siteIndex][ii];
      ofstreamExportFwdProb_ << ((ii < (probabilities[siteIndex].size() - 1)) ? "\t" : "\n");

    }

    siteIndex++;

  }

  ofstreamExportFwdProb_.close();

}

