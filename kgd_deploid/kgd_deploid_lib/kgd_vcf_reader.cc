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

#include "kgd_global.h"
#include <cassert>       // assert
#include <iostream>      // std::cout
#include "kgd_vcf_reader.h"


namespace kgd = kellerberrin::deploid;


/*! Initialize vcf file, search for the end of the vcf header.
 *  Extract the first block of data ( "buffer_length" lines ) into buff
 */
kgd::VcfReader::VcfReader(std::string fileName) {

  /*! Initialize by read in the vcf header file */
  init(fileName);
  readHeader();
  readVariants();
  getChromList();
  getIndexOfChromStarts();

  assert(doneGetIndexOfChromStarts_);

}


void kgd::VcfReader::checkFileCompressed() {

  FILE *f = NULL;

  f = fopen(fileName_.c_str(), "rb");

  if (f == NULL) {

    throw InvalidInputFile(fileName_);

  }

  unsigned char magic[2];

  fread((void *) magic, 1, 2, f);
  setIsCompressed((int(magic[0]) == 0x1f) && (int(magic[1]) == 0x8b));
  fclose(f);

}

void kgd::VcfReader::init(std::string fileName) {

  /*! Initialize other VcfReader class members
   */
  fileName_ = fileName;

  checkFileCompressed();

  if (isCompressed()) {

    inFileGz.open(fileName_.c_str(), std::ios::in);

  } else {

    inFile.open(fileName_.c_str(), std::ios::in);

  }

}


void kgd::VcfReader::finalize() {

  for (size_t i = 0; i < variants.size(); i++) {

    refCount.push_back((double) variants[i].ref);
    altCount.push_back((double) variants[i].alt);

  }

  if (isCompressed()) {

    inFileGz.close();

  } else {

    inFile.close();
  }

}


void kgd::VcfReader::readHeader() {

  if (isCompressed()) {

    if (!inFileGz.good()) {

      throw InvalidInputFile(fileName_);

    }

  } else {

    if (!inFile.good()) {

      throw InvalidInputFile(fileName_);

    }

  }

  if (isCompressed()) {

    getline(inFileGz, tmpLine_);

  } else {

    getline(inFile, tmpLine_);

  }

  while (tmpLine_.size() > 0) {

    if (tmpLine_[0] == '#') {

      if (tmpLine_[1] == '#') {

        headerLines.push_back(tmpLine_);

        if (isCompressed()) {

          getline(inFileGz, tmpLine_);

        } else {

          getline(inFile, tmpLine_);

        }

      } else {

        checkFeilds();
        break; //end of the header

      }

    } else {

      checkFeilds();

    }

  }


  dout << " There are " << headerLines.size() << " lines in the header." << std::endl;

}


void kgd::VcfReader::checkFeilds() {

  size_t field_start = 0;
  size_t field_end = 0;
  size_t field_index = 0;

  while (field_end < tmpLine_.size()) {

    field_end = std::min(tmpLine_.find('\t', field_start), tmpLine_.find('\n', field_start));
    tmpStr_ = tmpLine_.substr(field_start, field_end - field_start);
    std::string correctFieldValue;

    switch (field_index) {
      case 0:
        correctFieldValue = "#CHROM";
        break;
      case 1:
        correctFieldValue = "POS";
        break;
      case 2:
        correctFieldValue = "ID";
        break;
      case 3:
        correctFieldValue = "REF";
        break;
      case 4:
        correctFieldValue = "ALT";
        break;
      case 5:
        correctFieldValue = "QUAL";
        break;
      case 6:
        correctFieldValue = "FILTER";
        break;
      case 7:
        correctFieldValue = "INFO";
        break;
      case 8:
        correctFieldValue = "FORMAT";
        break;
      case 9:
        sampleName = tmpStr_;
        break;
    }

    if (tmpStr_ != correctFieldValue && field_index < 9) {
      throw VcfInvalidHeaderFieldNames(correctFieldValue, tmpStr_);
    }

    if (field_index == 9) {
      break;
    }

    field_start = field_end + 1;
    field_index++;

  } // End of while loop: field_end < line.size()

  assert(field_index == 9);
  assert(printSampleName());
}


void kgd::VcfReader::readVariants() {

  if (isCompressed()) {

    getline(inFileGz, tmpLine_);

  } else {

    getline(inFile, tmpLine_);

  }

  while (inFile.good() && tmpLine_.size() > 0) {

    VariantLine newVariant(tmpLine_);
    // check variantLine quality
    variants.push_back(newVariant);

    if (isCompressed()) {

      getline(inFileGz, tmpLine_);

    } else {

      getline(inFile, tmpLine_);

    }

  }

}


void kgd::VcfReader::getChromList() {

  chrom_.clear();
  position_.clear();

  assert (chrom_.size() == (size_t) 0);
  assert (position_.size() == (size_t) 0);

  std::string previousChrom("");
  std::vector<int> positionOfChrom_;

  for (size_t i = 0; i < variants.size(); i++) {

    if (previousChrom != variants[i].chromStr && previousChrom.size() > (size_t) 0) {

      chrom_.push_back(previousChrom);
      position_.push_back(positionOfChrom_);
      positionOfChrom_.clear();

    }

    positionOfChrom_.push_back(std::stoi(variants[i].posStr.c_str(), NULL));
    previousChrom = variants[i].chromStr;

  }

  chrom_.push_back(previousChrom);
  position_.push_back(positionOfChrom_);
  assert (position_.size() == chrom_.size());

}


void kgd::VcfReader::removeMarkers() {

  assert (keptVariants.size() == (size_t) 0);

  for (auto const &value: indexOfContentToBeKept) {

    keptVariants.push_back(variants[value]);

  }

  variants.clear();
  variants = keptVariants;
  keptVariants.clear();
  nLoci_ = variants.size();
  dout << " Vcf number of loci kept = " << nLoci_ << std::endl;

}


kgd::VariantLine::VariantLine(std::string tmpLine) {
  init(tmpLine);

  while (fieldEnd_ < tmpLine_.size()) {

    fieldEnd_ = std::min(tmpLine_.find('\t', feildStart_), tmpLine_.find('\n', feildStart_));
    tmpStr_ = tmpLine_.substr(feildStart_, fieldEnd_ - feildStart_);

    switch (fieldIndex_) {
      case 0:
        extract_field_CHROM();
        break;
      case 1:
        extract_field_POS();
        break;
      case 2:
        extract_field_ID();
        break;
      case 3:
        extract_field_REF();
        break;
      case 4:
        extract_field_ALT();
        break;
      case 5:
        extract_field_QUAL();
        break;
      case 6:
        extract_field_FILTER();
        break;
      case 7:
        extract_field_INFO();
        break;
      case 8:
        extract_field_FORMAT();
        break;
      case 9:
        extract_field_VARIANT();
        break;
    }

    feildStart_ = fieldEnd_ + 1;
    fieldIndex_++;

  }

}


void kgd::VariantLine::init(std::string tmpLine) {

  tmpLine_ = tmpLine;
  feildStart_ = 0;
  fieldEnd_ = 0;
  fieldIndex_ = 0;
  adFieldIndex_ = -1;

}


void kgd::VariantLine::extract_field_CHROM() {

  chromStr = tmpStr_;

}


void kgd::VariantLine::extract_field_POS() {

  posStr = tmpStr_;

}


void kgd::VariantLine::extract_field_ID() {

  idStr = tmpStr_;

}


void kgd::VariantLine::extract_field_REF() {

  refStr = tmpStr_;

}


void kgd::VariantLine::extract_field_ALT() {

  altStr = tmpStr_;

}


void kgd::VariantLine::extract_field_QUAL() { // Check for PASS

  qualStr = tmpStr_;

}


void kgd::VariantLine::extract_field_FILTER() {

  filterStr = tmpStr_;

}

void kgd::VariantLine::extract_field_INFO() {

  infoStr = tmpStr_;

}


void kgd::VariantLine::extract_field_FORMAT() {

  formatStr = tmpStr_;

  size_t field_start = 0;
  size_t field_end = 0;
  size_t field_index = 0;

  while (field_end < formatStr.size()) {

    field_end = std::min(formatStr.find(':', field_start), formatStr.find('\n', field_start));

    if ("AD" == formatStr.substr(field_start, field_end - field_start)) {

      adFieldIndex_ = field_index;
      break;

    }
    field_start = field_end + 1;
    field_index++;

  }
  if (adFieldIndex_ == -1) {

    throw VcfCoverageFieldNotFound(tmpStr_);

  }

  assert (adFieldIndex_ > -1);

}


void kgd::VariantLine::extract_field_VARIANT() {

  size_t field_start = 0;
  size_t field_end = 0;
  int field_index = 0;

  while (field_end < tmpStr_.size()) {

    field_end = std::min(tmpStr_.find(':', field_start), tmpStr_.find('\n', field_start));

    if (field_index == adFieldIndex_) {

      std::string adStr = tmpStr_.substr(field_start, field_end - field_start);
      size_t commaIndex = adStr.find(',', 0);
      ref = std::stoi(adStr.substr(0, commaIndex));
      alt = std::stoi(adStr.substr(commaIndex + 1, adStr.size()));
      break;

    }

    field_start = field_end + 1;
    field_index++;

  }

}


