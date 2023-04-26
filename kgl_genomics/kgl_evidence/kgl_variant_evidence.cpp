//
// Created by kellerberrin on 22/08/18.
//


#include "kgl_variant_evidence.h"
#include "kgl_data_file_type.h"

namespace kgl = kellerberrin::genome;


std::string kgl::VariantEvidence::output(char delimiter, VariantOutputIndex) const {

  std::stringstream ss;

  ss << delimiter << "VCFSource:" << DataDB::dataSource(dataSource());

  ss << delimiter   << "VCFRecord:" << vcfRecordCount();

  ss << delimiter << "Pass Filter:" << (passFilter() ? "true" : "false");

  ss << delimiter << "Alt Variants:" << altVariantCount();

  ss << delimiter << "Alt Variant Index:" << altVariantIndex();

  if (formatData()) {

    ss << formatData().value()->output(delimiter);

  } else {

    ss << delimiter << "No Format Data";

  }

  if (infoData()) {

    ss << delimiter << "Info Data Available";

  } else {

    ss << delimiter << "No Info Data";

  }


  return ss.str();

}


std::string kgl::FormatData::output(char delimiter) const {

  std::stringstream ss;

  ss << delimiter << "RefCount:" << refCount();
  ss << delimiter << "AltCount:" << altCount();
  ss << delimiter << "DPCount:" << DPCount();
  ss << delimiter << "GQProbWrongVariant:" << GQProbWrongVariant();
  ss << delimiter << "Quality:" << Quality();

  return ss.str();

}

