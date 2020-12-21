//
// Created by kellerberrin on 22/08/18.
//


#include "kgl_variant_evidence.h"

namespace kgl = kellerberrin::genome;


std::string kgl::VariantEvidence::output(char delimiter, VariantOutputIndex) const {

  std::stringstream ss;

  ss << delimiter   << "VCFRecord:" << vcfRecordCount();

  if (formatData()) {

    ss << formatData().value()->output(delimiter);

  }

  return ss.str();

}


std::string kgl::FormatData::output(char delimiter) const {

  std::stringstream ss;

  ss << delimiter << "RefCount:" << refCount();
  ss << delimiter << "AltCount:" << altCount();
  ss << delimiter << "DPCount:" << DPCount();
  ss << delimiter << "GQValue:" << GQValue();
  ss << delimiter << "Quality:" << Quality();

  return ss.str();

}

