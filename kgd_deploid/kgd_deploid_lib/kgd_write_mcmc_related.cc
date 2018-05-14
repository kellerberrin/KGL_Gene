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

#include <kgd_deploid_app.h>
#include "kgd_deploid_io.h"
#include "kgd_mcmc.h"


namespace kgl = kellerberrin::genome;
namespace kgd = kellerberrin::deploid;


void kgd::DEploidIO::writeMcmcRelated(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {


  writeProp(mcmcSample, useIBD);
  writeLLK(mcmcSample, useIBD);
  writeHap(mcmcSample, useIBD);

  
  if (useIBD == false) {

    writeVcf(mcmcSample);
    siteOfTwoSwitchOne_ = mcmcSample->siteOfTwoSwitchOne;
    siteOfTwoMissCopyOne_ = mcmcSample->siteOfTwoMissCopyOne;
    siteOfTwoSwitchTwo_ = mcmcSample->siteOfTwoSwitchTwo;
    siteOfTwoMissCopyTwo_ = mcmcSample->siteOfTwoMissCopyTwo;
    siteOfOneSwitchOne_ = mcmcSample->siteOfOneSwitchOne;
    siteOfOneMissCopyOne_ = mcmcSample->siteOfOneMissCopyOne;

    finalSiteOfTwoSwitchOne_ = mcmcSample->currentsiteOfTwoSwitchOne;
    finalSiteOfTwoMissCopyOne_ = mcmcSample->currentsiteOfTwoMissCopyOne;
    finalSiteOfTwoSwitchTwo_ = mcmcSample->currentsiteOfTwoSwitchTwo;
    finalSiteOfTwoMissCopyTwo_ = mcmcSample->currentsiteOfTwoMissCopyTwo;
    finalSiteOfOneSwitchOne_ = mcmcSample->currentsiteOfOneSwitchOne;
    finalSiteOfOneMissCopyOne_ = mcmcSample->currentsiteOfOneMissCopyOne;

    //writeEventCount( );
  } else {
  
    //IBD_path_change_at_ = mcmcSample->IBD_path_change_at_;
    //finalIBDpathChangeAt = mcmcSample->current_IBD_path_change_at_;

  }

}


void kgd::DEploidIO::writeProp(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {

  if (useIBD) {

    ofstreamExportTmp_.open(strIbdExportProp_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  } else {

    ofstreamExportTmp_.open(strExportProp_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  }

  for (size_t i = 0; i < mcmcSample->proportion.size(); i++) {

    for (size_t ii = 0; ii < mcmcSample->proportion[i].size(); ii++) {

      ofstreamExportTmp_ << std::setw(10) << mcmcSample->proportion[i][ii];
      ofstreamExportTmp_ << ((ii < (mcmcSample->proportion[i].size() - 1)) ? "\t" : "\n");

    }

  }

  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeLLK(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {

  if (useIBD) {

    ofstreamExportTmp_.open(strIbdExportLLK_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  } else {

    ofstreamExportTmp_.open(strExportLLK_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  }

  for (size_t i = 0; i < mcmcSample->sumLLKs.size(); i++) {

    ofstreamExportTmp_ << mcmcSample->moves[i] << "\t" << mcmcSample->sumLLKs[i] << std::endl;

  }

  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeHap(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {

  if (useIBD) {

    ofstreamExportTmp_.open(strIbdExportHap_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  } else {

    ofstreamExportTmp_.open(strExportHap_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  }
  // HEADER
  ofstreamExportTmp_ << "CHROM" << "\t" << "POS" << "\t";;

  for (size_t ii = 0; ii < kStrain_; ii++) {

    ofstreamExportTmp_ << "h" << (ii + 1);
    ofstreamExportTmp_ << ((ii < (kStrain_ - 1)) ? "\t" : "\n");

  }

  size_t siteIndex = 0;

  for (size_t chromI = 0; chromI < chrom_.size(); chromI++) {

    for (size_t posI = 0; posI < position_[chromI].size(); posI++) {

      ofstreamExportTmp_ << chrom_[chromI] << "\t" << (int) position_[chromI][posI] << "\t";

      for (size_t ii = 0; ii < mcmcSample->hap[siteIndex].size(); ii++) {

        ofstreamExportTmp_ << mcmcSample->hap[siteIndex][ii];
        ofstreamExportTmp_ << ((ii < (mcmcSample->hap[siteIndex].size() - 1)) ? "\t" : "\n");

      }

      siteIndex++;

    }

  }

  assert (siteIndex == mcmcSample->hap.size());
  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeVcf(std::shared_ptr<McmcSample> mcmcSample) {

  if (!doExportVcf()) return;

  ogzstream ogstreamExport;
  std::ostream *writeTo;

  if (compressVcf()) {

    ogstreamExport.open(strExportVcf_.c_str(), std::ios::out);
    writeTo = &ogstreamExport;

  } else {

    ofstreamExportTmp_.open(strExportVcf_.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    writeTo = &ofstreamExportTmp_;

  }

  // VCF HEADER
  if (useVcf()) {

    for (auto const &headerLine: vcfReaderPtr_->headerLines) {

      (*writeTo) << headerLine << std::endl;

    }

  } else {

    (*writeTo) << "##fileformat=VCFv4.2" << std::endl;

  }
  // DEploid call

  (*writeTo) << "##DEploid call: kgd_deploid ";
  (*writeTo) << kgl::ExecEnv::commandLine() << std::endl;

  // Include proportions
  for (size_t ii = 0; ii < kStrain_; ii++) {

    (*writeTo) << "##Proportion of strain "
               << (useVcf() ? vcfReaderPtr_->sampleName : "h")
               << "." << (ii + 1)
               << "=" << mcmcSample->proportion.back()[ii] << std::endl;

  }

  // HEADER
  (*writeTo) << "#CHROM" << "\t"
             << "POS" << "\t"
             << "ID" << "\t"
             << "REF" << "\t"
             << "ALT" << "\t"
             << "QUAL" << "\t"
             << "FILTER" << "\t"
             << "INFO" << "\t"
             << "FORMAT" << "\t";

  for (size_t ii = 0; ii < kStrain_; ii++) {

    (*writeTo) << (useVcf() ? vcfReaderPtr_->sampleName : "h")
               << "." << (ii + 1);
    (*writeTo) << ((ii < (kStrain_ - 1)) ? "\t" : "\n");

  }

  size_t siteIndex = 0;

  for (size_t chromI = 0; chromI < chrom_.size(); chromI++) {

    for (size_t posI = 0; posI < position_[chromI].size(); posI++) {

      if (useVcf()) {

        (*writeTo) << vcfReaderPtr_->variants[siteIndex].chromStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].posStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].idStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].refStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].altStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].qualStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].filterStr << "\t"
                   << vcfReaderPtr_->variants[siteIndex].infoStr << "\t"
                   << "GT" << "\t";
      } else {

        (*writeTo) << chrom_[chromI] << "\t"
                   << (int) position_[chromI][posI] << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "GT" << "\t";
      }

      for (size_t ii = 0; ii < mcmcSample->hap[siteIndex].size(); ii++) {

        (*writeTo) << mcmcSample->hap[siteIndex][ii];
        (*writeTo) << ((ii < (mcmcSample->hap[siteIndex].size() - 1)) ? "\t" : "\n");

      }

      siteIndex++;

    }

  }

  assert (siteIndex == mcmcSample->hap.size());

  if (compressVcf()) {

    ogstreamExport.close();

  } else {

    ofstreamExportTmp_.close();

  }

}

