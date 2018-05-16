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


#ifndef KGD_VCF_H
#define KGD_VCF_H



#include <string>  /* string */
#include <vector>  /* vector */
#include <fstream>
#include <stdlib.h>     /* strtol, strtod */
#include "kgd_exceptions.h"
#include "kgd_variant_index.h"
#include "gzstream.h"



namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


  
struct InvalidVcf : public InvalidInput {

  InvalidVcf(std::string str) : InvalidInput(str) {}

  virtual ~InvalidVcf() throw() {}
  //virtual const char* what () const noexcept {
  //return throwMsg.c_str();
  //}
};


struct VcfInvalidHeaderFieldNames : public InvalidVcf {

  VcfInvalidHeaderFieldNames(std::string str1, std::string str2) : InvalidVcf(str1) {

    this->reason = " VCF field header expects: ";
    throwMsg = this->reason + this->src + ", " + str2 + " was found!";

  }

  ~VcfInvalidHeaderFieldNames() throw() {}
};


struct VcfInvalidVariantEntry : public InvalidVcf {

  VcfInvalidVariantEntry(std::string str) : InvalidVcf(str) {}

  virtual ~VcfInvalidVariantEntry() throw() {}
  //virtual const char* what () const noexcept {
  //return throwMsg.c_str();
  //}
};


struct VcfCoverageFieldNotFound : public VcfInvalidVariantEntry {

  VcfCoverageFieldNotFound(std::string str) : VcfInvalidVariantEntry(str) {

    this->reason = "Coverage field AD was not found in the FORMAT, found: ";
    throwMsg = this->reason + this->src;

  }

  ~VcfCoverageFieldNotFound() throw() {}

};

// This requires more thinking, check for data type?
// More informative exceptions for vcf related errors


class VariantLine {


public:


  VariantLine(std::string tmpLine);
  ~VariantLine() = default;

  // Access functions.
  double getFloatRef() const { return static_cast<double>(ref_); }
  double getFloatAlt() const { return static_cast<double>(alt_); }
  const std::string& getChromStr() const { return chromStr_; }
  const std::string& getPosStr() const { return posStr_; }
  const std::string& getIdStr() const { return idStr_; }
  const std::string& getRefStr() const { return refStr_; }
  const std::string& getAltStr() const { return altStr_; }
  const std::string& getQualStr() const { return qualStr_; }
  const std::string& getFilterStr() const { return filterStr_; }
  const std::string& getInfoStr() const { return infoStr_; }
  const std::string& getFormatStr() const { return formatStr_; }

private:

  std::string tmpLine_;
  std::string tmpStr_;

  size_t feildStart_;
  size_t fieldEnd_;
  size_t fieldIndex_;

  std::string chromStr_;
  std::string posStr_;
  std::string idStr_;
  std::string refStr_;
  std::string altStr_;
  std::string qualStr_;
  std::string filterStr_;
  std::string infoStr_;
  std::string formatStr_;

  int adFieldIndex_;
  int ref_;
  int alt_;

  void init(std::string tmpLine);

  void extract_field_CHROM();
  void extract_field_POS();
  void extract_field_ID();
  void extract_field_REF();
  void extract_field_ALT();
  void extract_field_QUAL();
  void extract_field_FILTER();
  void extract_field_INFO();
  void extract_field_FORMAT();
  void extract_field_VARIANT();


};


/*! \brief VCF file reader @ingroup group_data */
class VcfReader : public VariantIndex {

#ifdef UNITTEST
  friend class TestVCF;
#endif


public:

  // Constructors and Destructors
  VcfReader(std::string fileName); // parse in exclude sites
  ~VcfReader() = default;

  // Access functions.
  const std::vector<std::string>& getHeaderLines() const { return headerLines_; }
  const std::string& getSampleName() const { return sampleName_; }
  const VariantLine& getVariantLine(size_t index) const { return variants_[index]; }
  const std::vector<double>& getRefCount() const { return refCount_; }
  const std::vector<double>& getAltCount() const { return altCount_; }

  // Methods
  void finalize();

private:

  std::vector<VariantLine> variants_;
  std::vector<VariantLine> keptVariants_;
  std::vector<double> refCount_;
  std::vector<double> altCount_;
  std::vector<std::string> headerLines_;
  std::string fileName_;

  std::string sampleName_;
  std::string tmpLine_;
  std::string tmpStr_;

  std::ifstream inFile_;
  igzstream inFileGz;
  bool isCompressed_;

  // Getters ands Setters.
  bool isCompressed() const { return this->isCompressed_; }
  void setIsCompressed(const bool compressed) { this->isCompressed_ = compressed; }
  void checkFileCompressed();

  // Methods
  void init(std::string fileName);
  void readVariants();
  void readHeader();
  void checkFeilds();
  void getChromList();
  void removeMarkers();

  // Debug tools
  bool printSampleName();


};



}   // organization level namespace
}   // project level namespace



#endif
