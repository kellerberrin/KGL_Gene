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

#include "kgd_dEploidIO.h"
#include "kgd_mcmc.h"
#include "kgd_deploid_app.h"


namespace kgd = kellerberrin::deploid;


void kgd::DEploidIO::wrapUp() {

  writeRecombProb(panel);

  // Get End time before writing the log
  getTime(false);

  writeLog(&std::cout);

  ofstreamExportTmp_.open(strExportLog_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  writeLog(&ofstreamExportTmp_);

  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeRecombProb(std::shared_ptr<Panel> panel) {

  if (!doExportRecombProb()) return;

  if (panel) {

    ofstreamExportTmp_.open(strExportRecombProb_.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    ofstreamExportTmp_ << "p.recomb" << "\t"
                      << "p.each" << "\t"
                      << "p.no.recomb" << "\t"
                      << "p.rec.rec" << "\t"
                      << "p.rec.norec" << "\t"
                      << "p.norec.norec" << "\n";

    for (size_t i = 0; i < panel->pRec_.size(); i++) {

      ofstreamExportTmp_ << panel->pRec_[i] << "\t"
                        << panel->pRecEachHap_[i] << "\t"
                        << panel->pNoRec_[i] << "\t"
                        << panel->pRecRec_[i] << "\t"
                        << panel->pRecNoRec_[i] << "\t"
                        << panel->pNoRecNoRec_[i] << "\n";
    }

    ofstreamExportTmp_.close();

  }

}


void kgd::DEploidIO::writeLog(std::ostream *writeTo) {

  size_t nHash = 30 + std::string(DeploidExecEnv::VERSION).size();

  for (size_t i = 0; i < nHash; i++) {

    (*writeTo) << "#";

  }

  (*writeTo) << "\n";
  (*writeTo) << "#        kgd_deploid " << std::setw(10) << DeploidExecEnv::VERSION << " log        #\n";

  for (size_t i = 0; i < nHash; i++) {

    (*writeTo) << "#";

  }

  (*writeTo) << "\n";
  (*writeTo) << "Program was compiled on: " << compileTime_ << std::endl;
  (*writeTo) << "kgd_deploid version: " << dEploidGitVersion_ << std::endl;
  (*writeTo) << "\n";
  (*writeTo) << "Input data: \n";

  if (panelFileName_.size() > 0) {

    (*writeTo) << std::setw(12) << "Panel: " << panelFileName_ << "\n";

  }

  (*writeTo) << std::setw(12) << "PLAF: " << plafFileName_ << "\n";

  if (useVcf()) (*writeTo) << std::setw(12) << "VCF: " << vcfFileName_ << "\n";

  if (refFileName_.size() > 0) (*writeTo) << std::setw(12) << "REF count: " << refFileName_ << "\n";

  if (altFileName_.size() > 0) (*writeTo) << std::setw(12) << "ALT count: " << altFileName_ << "\n";

  if (excludeSites()) { (*writeTo) << std::setw(12) << "Exclude: " << excludeFileName_ << "\n"; }

  (*writeTo) << "\n";

  if ((doLsPainting() == false) & (doIbdPainting() == false)) {

    (*writeTo) << "MCMC parameters: " << "\n";
    (*writeTo) << std::setw(19) << " MCMC burn: " << mcmcBurn_ << "\n";
    (*writeTo) << std::setw(19) << " MCMC sample: " << nMcmcSample_ << "\n";
    (*writeTo) << std::setw(19) << " MCMC sample rate: " << mcmcMachineryRate_ << "\n";
    (*writeTo) << std::setw(19) << " Random seed: " << randomSeed() << "\n";

    if (useIBD()) {

      (*writeTo) << std::setw(19) << "  IBD Method used: YES" << "\n";

    }

    (*writeTo) << std::setw(19) << " Update Prop: " << (doUpdateProp() ? "YES" : "NO") << "\n";
    (*writeTo) << std::setw(19) << " Update Single: " << (doUpdateSingle() ? "YES" : "NO") << "\n";
    (*writeTo) << std::setw(19) << " Update Pair: " << (doUpdatePair() ? "YES" : "NO") << "\n";
    (*writeTo) << "\n";

  }

  (*writeTo) << "Other parameters:" << "\n";

  if (forbidCopyFromSame_) {

    (*writeTo) << " Update pair haplotypes move forbid copying from the same strain!!! \n";

  }

  (*writeTo) << std::setw(20) << " Miss copy prob: " << missCopyProb_ << "\n";
  (*writeTo) << std::setw(20) << " Avrg Cent Morgan: " << averageCentimorganDistance_ << "\n";
  (*writeTo) << std::setw(20) << " G: " << parameterG() << "\n";

  if (useIBD()) {

    (*writeTo) << std::setw(20) << " IBD sigma: " << ibdSigma() << "\n";

  } else {

    (*writeTo) << std::setw(20) << " sigma: " << parameterSigma() << "\n";

  }

  (*writeTo) << std::setw(20) << " ScalingFactor: " << scalingFactor() << "\n";

  if (initialPropWasGiven()) {

    (*writeTo) << std::setw(20) << " Initial prob: ";

    for (size_t i = 0; i < initialProp_.size(); i++) {

      (*writeTo) << initialProp_[i] << ((i != (kStrain_ - 1)) ? " " : "\n");

    }

  }

  (*writeTo) << "\n";

  if ((doLsPainting() == false) & (doIbdPainting() == false) & (doComputeLLK() == false)) {

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

  if (doComputeLLK()) {

    (*writeTo) << "Input likelihood: " << llkFromInitialHap_;
    (*writeTo) << "\n";

  } else {

    (*writeTo) << "Output saved to:\n";

    if (doLsPainting()) {

      for (size_t i = 0; i < kStrain(); i++) {

        (*writeTo) << "Posterior probability of strain " << i << ": " << strExportSingleFwdProbPrefix_ << i << std::endl;

      }

    } else if (doIbdPainting()) {

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

      if (doExportVcf()) { (*writeTo) << std::setw(14) << "Vcf: " << strExportVcf_ << "\n"; }

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
                    << "IBDpathChangeAt" << "\t"
                    << "finalIBDpathChangeAt" << "\t"

                    << "siteOfTwoSwitchOne" << "\t"
                    << "finalSiteOfTwoSwitchOne" << "\t"

                    << "siteOfTwoMissCopyOne" << "\t"
                    << "finalSiteOfTwoMissCopyOne" << "\t"

                    << "siteOfTwoSwitchTwo" << "\t"
                    << "finalSiteOfTwoSwitchTwo" << "\t"

                    << "siteOfTwoMissCopyTwo" << "\t"
                    << "finalSiteOfTwoMissCopyTwo" << "\t"

                    << "siteOfOneSwitchOne" << "\t"
                    << "finalSiteOfOneSwitchOne" << "\t"

                    << "siteOfOneMissCopyOne" << "\t"
                    << "finalSiteOfOneMissCopyOne" << std::endl;

  size_t siteIndex = 0;

  for (size_t chromI = 0; chromI < chrom_.size(); chromI++) {

    for (size_t posI = 0; posI < position_[chromI].size(); posI++) {

      ofstreamExportTmp_ << chrom_[chromI] << "\t"
                        << (int) position_[chromI][posI] << "\t"

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


void kgd::DEploidIO::writeIBDpostProb(std::vector<std::vector<double> > &reshapedProbs, std::vector<std::string> header) {

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

  for (size_t chromIndex = 0; chromIndex < position_.size(); chromIndex++) {

    for (size_t posI = 0; posI < position_[chromIndex].size(); posI++) {

      (*writeTo) << chrom_[chromIndex] << "\t" << (int) position_[chromIndex][posI] << "\t";

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
