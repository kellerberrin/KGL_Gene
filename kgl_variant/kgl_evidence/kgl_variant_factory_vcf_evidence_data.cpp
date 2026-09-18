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

  return {InfoEvidenceSubscriber::NotImplemented, InfoEvidenceExtern::Boolean, InfoEvidenceIntern::NotImplemented };

}


