//
// Created by kellerberrin on 23/5/20.
//

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_data.h"



namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper class for INFO types.


// Returns true if a number > 1. Used in the type definition lambdas.
bool kgl::InfoTypeLookup::isVectorType(const std::string& type) {

  try {

    return std::stol(type) > 1;

  }
  catch(...) {

    return false;

  }

}


kgl::InfoEvidenceType kgl::InfoTypeLookup::evidenceType(const VCFInfoRecord &vcf_info_item) {

  for (auto const& evidence_type : type_definitions_) {

    if (evidence_type.second(vcf_info_item.type, vcf_info_item.number)) {

      return evidence_type.first;

    }

  }

  ExecEnv::log().warn("InfoTypeLookup::evidenceType, Info ID: {} Unable to find data type combination for Number: {}, Type: {}",
                      vcf_info_item.ID, vcf_info_item.number, vcf_info_item.type);

  return {InfoEvidenceSubscriber::NotImplemented, InfoEvidenceExtern::NotImplemented, InfoEvidenceIntern::NotImplemented };

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

bool kgl::InfoEvidenceType::fixedDataType() const {

  switch (InternalInfoType()) {

    case InfoEvidenceIntern::intern_char: // Size is fixed (boolean variables) and known at Info data subscription time.
    case InfoEvidenceIntern::intern_integer: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::intern_float: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::intern_integer_array: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::intern_float_array: // Size is fixed and known at Info data subscription time.
    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
      return true;

    case InfoEvidenceIntern::intern_string: // Size varies between records.
    case InfoEvidenceIntern::intern_string_array:  // Size varies between records.
    case InfoEvidenceIntern::intern_unity_integer_array:   // Size varies between records.
    case InfoEvidenceIntern::intern_unity_float_array:   // Size varies between records.
    case InfoEvidenceIntern::intern_unity_string_array:    // Size varies between records.
      return false;

    default:
      return false;
  }

}



