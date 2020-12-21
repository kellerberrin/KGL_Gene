
#include <iomanip>      // std::setw
#include <kgd_deconvolv_app.h>
#include "kgd_deploid_io.h"
#include "kgd_mcmc_sample.h"


namespace kgd = kellerberrin::deconvolv;


void kgd::DEploidIO::writeMcmcRelated(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {


  writeProp(mcmcSample, useIBD);
  writeLLK(mcmcSample, useIBD);
  writeHap(mcmcSample, useIBD);

  if (not useIBD) {

    writeVcf(mcmcSample);

  }

}


void kgd::DEploidIO::writeProp(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {

  std::string filename;

  if (useIBD) {

    filename = strIbdExportProp_;

  } else {

    filename = strExportProp_;

  }

  filename = "/home/kellerberrin/test.txt";

  std::ofstream fExport(filename);

  if (fExport.good()) {

    ExecEnv::log().info("Write strain proportions to file: {}", filename);

  } else {

    ExecEnv::log().warn("Unable to write to file: {}", filename);
    return;

  }

  size_t count = 0;
  for (size_t i = 0; i < mcmcSample->getProportion().size(); i++) {

    for (size_t ii = 0; ii < mcmcSample->getProportionIndex(i).size(); ii++) {

      fExport << std::setw(10) << mcmcSample->getProportionIndex(i, ii);
      fExport << ((ii < (mcmcSample->getProportionIndex(i).size() - 1)) ? "\t" : "\n");
      ++count;

    }

  }

  ExecEnv::log().info("{} strain proportions written", count);

}


void kgd::DEploidIO::writeLLK(std::shared_ptr<McmcSample> mcmcSample, bool useIBD) {

  if (useIBD) {

    ofstreamExportTmp_.open(strIbdExportLLK_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  } else {

    ofstreamExportTmp_.open(strExportLLK_.c_str(), std::ios::out | std::ios::app | std::ios::binary);

  }

  for (size_t i = 0; i < mcmcSample->getSumLLKs().size(); i++) {

    ofstreamExportTmp_ << mcmcSample->getMovesIndex(i) << "\t" << mcmcSample-> getSumLLKsIndex(i) << std::endl;

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

  for (size_t ii = 0; ii < kStrain(); ii++) {

    ofstreamExportTmp_ << "h" << (ii + 1);
    ofstreamExportTmp_ << ((ii < (kStrain() - 1)) ? "\t" : "\n");

  }

  size_t siteIndex = 0;

  for (size_t chromI = 0; chromI < getMixtureData().getChrom().size(); chromI++) {

    for (size_t posI = 0; posI < getMixtureData().getPosition()[chromI].size(); posI++) {

      ofstreamExportTmp_ << getMixtureData().getChrom()[chromI] << "\t" << static_cast<int>(getMixtureData().getPosition()[chromI][posI]) << "\t";

      for (size_t ii = 0; ii < mcmcSample->getHapIndex(siteIndex).size(); ii++) {

        ofstreamExportTmp_ << mcmcSample->getHapIndex(siteIndex, ii);
        ofstreamExportTmp_ << ((ii < (mcmcSample->getHapIndex(siteIndex).size() - 1)) ? "\t" : "\n");

      }

      siteIndex++;

    }

  }

  assert (siteIndex == mcmcSample->getHap().size());
  ofstreamExportTmp_.close();

}


void kgd::DEploidIO::writeVcf(std::shared_ptr<McmcSample> mcmcSample) {

  if (not getMixtureControl().doExportVcf()) return;

  if (not getMixtureControl().useVcf()) {

    ExecEnv::log().warn("Can only export VCF if VCF file specified");
    return;

  }

  ogzstream ogstreamExport;
  std::ostream *writeTo;

  if (getMixtureControl().compressVcf()) {

    ogstreamExport.open(strExportVcf_.c_str(), std::ios::out);
    writeTo = &ogstreamExport;

  } else {

    ofstreamExportTmp_.open(strExportVcf_.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    writeTo = &ofstreamExportTmp_;

  }

  VcfReader VCF_reader(vcfFileName_);

  // VCF HEADER
  if (getMixtureControl().useVcf()) {

    for (auto const &headerLine: VCF_reader.getHeaderLines()) {

      (*writeTo) << headerLine << std::endl;

    }

  } else {

    (*writeTo) << "##fileformat=VCFv4.2" << std::endl;

  }
  // DEploid call

  (*writeTo) << "##DEploid call: kgd_deconvolv ";
  (*writeTo) << ExecEnv::commandLine() << std::endl;

  // Include proportions
  for (size_t ii = 0; ii < kStrain(); ii++) {

    std::string sampleName = VCF_reader.getSampleName();

    (*writeTo) << "##Proportion of strain "
               << (getMixtureControl().useVcf() ? sampleName : "h")
               << "." << (ii + 1)
               << "=" << mcmcSample->getProportionIndex(ii).back() << std::endl;

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

  for (size_t ii = 0; ii < kStrain(); ii++) {

    std::string sampleName = VCF_reader.getSampleName();

    (*writeTo) << (getMixtureControl().useVcf() ? sampleName : "h")
               << "." << (ii + 1);
    (*writeTo) << ((ii < (kStrain() - 1)) ? "\t" : "\n");

  }

  size_t siteIndex = 0;

  for (size_t chromI = 0; chromI < getMixtureData().getChrom().size(); chromI++) {

    for (size_t posI = 0; posI < getMixtureData().getPosition()[chromI].size(); posI++) {

      if (getMixtureControl().useVcf()) {

        (*writeTo) << VCF_reader.getVariantLine(siteIndex).getChromStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getPosStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getIdStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getRefStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getAltStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getQualStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getFilterStr() << "\t"
                   << VCF_reader.getVariantLine(siteIndex).getInfoStr() << "\t"
                   << "GT" << "\t";
      } else {

        (*writeTo) << getMixtureData().getChrom()[chromI] << "\t"
                   << static_cast<int>(getMixtureData().getPosition()[chromI][posI]) << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "." << "\t"
                   << "GT" << "\t";
      }

      for (size_t ii = 0; ii < mcmcSample->getHapIndex(siteIndex).size(); ii++) {

        (*writeTo) << mcmcSample->getHapIndex(siteIndex, ii);
        (*writeTo) << ((ii < (mcmcSample->getHapIndex(siteIndex).size() - 1)) ? "\t" : "\n");

      }

      siteIndex++;

    }

  }

  assert (siteIndex == mcmcSample->getHap().size());

  if (getMixtureControl().compressVcf()) {

    ogstreamExport.close();

  } else {

    ofstreamExportTmp_.close();

  }

}

