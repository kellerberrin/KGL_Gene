//
// Created by kellerberrin on 16/5/20.
//

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoEvidenceIndex. There is only one of these held by all variants with INFO evidence fields.

bool kgl::InfoEvidenceHeader::setupEvidenceHeader(const VCFInfoRecordMap& vcf_info_map, std::shared_ptr<const InfoEvidenceHeader> self_ptr) {

  // Create a map of indexes.
  info_index_map_.clear();
  for (auto const& [ident , info_record] : vcf_info_map) {

    InfoEvidenceType type = InfoTypeCount::convertVCFType(info_record);
    InfoEvidenceIndex info_index(info_record, type, 0, self_ptr);
    auto result = info_index_map_.emplace(ident, info_index);
    if (not result.second) {

      ExecEnv::log().error("InfoEvidenceHeader::setupEvidenceHeader, could not add index for info field: {}", ident);

    }

  }

  return true;
}


std::optional<const kgl::InfoEvidenceIndex> kgl::InfoEvidenceHeader::getIndex(const std::string&) const {

  return std::nullopt;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper class for INFO types.


std::optional<kgl::InfoIndexMap> kgl::InfoTypeCount::generateEvidenceIndex(const VCFInfoRecordMap& vcf_info_map) {


  return std::nullopt;
}


void  kgl::InfoTypeCount::reset() {

  float_count_ = 0;
  integer_count_ = 0;
  string_count_ = 0;
  float_array_count_ = 0;
  integer_array_count_ = 0;
  string_array_count_ = 0;
  boolean_count_ = 0;

}

void kgl::InfoTypeCount::incrementCount(InfoEvidenceType type) {

  switch(type) {

    case InfoEvidenceType::Float:
      ++float_count_;
      break;

    case InfoEvidenceType::Integer:
      ++integer_count_;
      break;

    case InfoEvidenceType::String:
      ++string_count_;
      break;

    case InfoEvidenceType::FloatArray:
      ++float_array_count_;
      break;

    case InfoEvidenceType::IntegerArray:
      ++integer_array_count_;
      break;

    case InfoEvidenceType::StringArray:
      ++string_array_count_;
      break;

    case InfoEvidenceType::Boolean:
      ++boolean_count_;
      break;

    case InfoEvidenceType::FloatAllele:
      ++float_allele_count_;
      break;

    case InfoEvidenceType::IntegerAllele:
      ++integer_allele_count_;
      break;

    case InfoEvidenceType::StringAllele:
      ++string_allele_count_;
      break;


    case InfoEvidenceType::NotImplemented:
    default:
      break;

  }

}


kgl::InfoEvidenceType kgl::InfoTypeCount::convertVCFType(const VCFInfoRecord& vcf_info_item) {

  try {

    if (vcf_info_item.type == INTEGER_) {

      if (vcf_info_item.number == SCALAR_) {

        return InfoEvidenceType::Integer;

      } else if (vcf_info_item.number == ALLELE_) {

        return InfoEvidenceType::IntegerAllele;

      } else if (vcf_info_item.number == INDETERMINATE_COUNT_) {

        return InfoEvidenceType::IntegerArray;

      } else if (std::stol(vcf_info_item.number) > 1) {

        return InfoEvidenceType::IntegerArray;

      } else {

        return InfoEvidenceType::NotImplemented;

      }


    } else if (vcf_info_item.type == FLOAT_) {

      if (vcf_info_item.number == SCALAR_) {

        return InfoEvidenceType::Float;

      } else if (vcf_info_item.number == ALLELE_) {

        return InfoEvidenceType::FloatAllele;

      } else if (vcf_info_item.number == INDETERMINATE_COUNT_) {

        return InfoEvidenceType::FloatArray;

      } else if (std::stol(vcf_info_item.number) > 1) {

        return InfoEvidenceType::FloatArray;

      } else {

        ExecEnv::log().warn(
        "InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
        vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
        return InfoEvidenceType::NotImplemented;

      }

    } else if (vcf_info_item.type == FLAG_) {

      if (vcf_info_item.number == FLAG_SCALAR_) {

        return InfoEvidenceType::Boolean;

      } else {

        ExecEnv::log().warn(
        "InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
        vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
        return InfoEvidenceType::NotImplemented;

      }

    } else if (vcf_info_item.type == CHAR_STRING_ or vcf_info_item.type == STRING_) {

      if (vcf_info_item.number == SCALAR_) {

        return InfoEvidenceType::String;

      } else if (vcf_info_item.number == ALLELE_) {

        return InfoEvidenceType::StringAllele;

      } else if (vcf_info_item.number == INDETERMINATE_COUNT_) {

        return InfoEvidenceType::StringArray;

      } else if (std::stol(vcf_info_item.number) > 1) {

        return InfoEvidenceType::StringArray;

      } else {

        ExecEnv::log().warn(
        "InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
        vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
        return InfoEvidenceType::NotImplemented;

      }

    } else {

      ExecEnv::log().warn(
      "InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
      vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
      return InfoEvidenceType::NotImplemented;

    }

  }
  catch(...) {

//    ExecEnv::log().error(
//    "std::stol Unknown Exception; InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
//    vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
    return InfoEvidenceType::NotImplemented;

  }

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant (data only).


kgl::InfoDataPair kgl::EvidenceFactory::createVariantEvidence(std::string&& info) {

  // If no Info fields have been requested, then just return std::nullopt
  if (evidence_map_.empty()) {

    return InfoDataPair(std::nullopt, std::nullopt);

  }

//  test1(info);
  test2(info);

  // Just return a single data block for now
  return InfoDataPair(std::make_unique<InfoDataBlock>(info_evidence_header_), std::nullopt);

}

void kgl::EvidenceFactory::test1(std::string& info) {


  // Parse the info line.
  VCFInfoParser info_parser_(std::move(info));

  for (auto const& subscribed_info : active_info_map_) {

    auto vcf_info_item = info_parser_.getInfoString(subscribed_info.first);

    if (vcf_info_item) {

      try {

        if (InfoTypeCount::convertVCFType(subscribed_info.second) == InfoEvidenceType::FloatAllele or
            InfoTypeCount::convertVCFType(subscribed_info.second) == InfoEvidenceType::Float) {

          std::stod(vcf_info_item.value());

        }

        if (InfoTypeCount::convertVCFType(subscribed_info.second) == InfoEvidenceType::IntegerAllele
            or InfoTypeCount::convertVCFType(subscribed_info.second) == InfoEvidenceType::IntegerAllele) {

          std::stol(vcf_info_item.value());

        }

        if (InfoTypeCount::convertVCFType(subscribed_info.second) == InfoEvidenceType::IntegerArray) {

          std::vector<std::string> int_array_vec = Utility::tokenizer(vcf_info_item.value(), ",");

          for (auto const& int_str : int_array_vec) {

            std::stol(int_str);

          }

        }

        if (InfoTypeCount::convertVCFType(subscribed_info.second) == InfoEvidenceType::FloatArray) {

          std::vector<std::string> float_array_vec = Utility::tokenizer(vcf_info_item.value(), ",");

          for (auto const& float_str : float_array_vec) {

            std::stod(float_str);

          }

        }

      }
      catch(std::out_of_range& e) {

        ExecEnv::log().error("Conversion Out of of Range Exception:{}, Field ID: {}, Desc: {}, Number: {}, Type: {}, Value: {}",
                             e.what(), subscribed_info.second.ID, subscribed_info.second.description,
                             subscribed_info.second.number, subscribed_info.second.type,
                             vcf_info_item.value());

      }
      catch(std::invalid_argument& e) {

        ExecEnv::log().error("Conversion Invalid Argument Exception:{}, Field ID: {}, Desc: {}, Number: {}, Type: {}, Value: {}",
                             e.what(), subscribed_info.second.ID, subscribed_info.second.description,
                             subscribed_info.second.number, subscribed_info.second.type,
                             vcf_info_item.value());

      }
      catch(std::exception& e) {

        ExecEnv::log().error("Conversion Unknown Exception, Field ID: {}, Desc: {}, Number: {}, Type: {}, Value: {}",
                             e.what(), subscribed_info.second.ID, subscribed_info.second.description,
                             subscribed_info.second.number, subscribed_info.second.type,
                             vcf_info_item.value());
      }


    } else {

      // Mark the bitfield as present.

    }

  }

}


void kgl::EvidenceFactory::test2(std::string& info) {


  // Parse the info line.
  VCFInfoParser info_parser_(std::move(info));

  for (auto const& subscribed_info : active_info_map_) {

    switch(InfoTypeCount::convertVCFType(subscribed_info.second)) {

      case InfoEvidenceType::FloatAllele:
      case InfoEvidenceType::Float:
        info_parser_.getInfoFloat(subscribed_info.first);
        break;

      case InfoEvidenceType::IntegerAllele:
      case InfoEvidenceType::Integer:
        break;

      case InfoEvidenceType::IntegerArray:
        info_parser_.getInfoIntegerArray(subscribed_info.first);
        break;

      case InfoEvidenceType::FloatArray:
        info_parser_.getInfoFloatArray(subscribed_info.first);
        break;

      default:
        break;

    }

  }

}


void kgl::EvidenceFactory::availableInfoFields(const VCFInfoRecordMap& vcf_info_map) {


  all_available_map_ = vcf_info_map;

  for (auto const& [ident, vcf_info_field] : vcf_info_map) {

    active_info_map_[ident] = vcf_info_field;

  }

/*
  for (auto const& subscribed_item : evidence_map_) {

    auto result = vcf_info_map.find(subscribed_item);

    if (result == vcf_info_map.end()) {

      ExecEnv::log().warn("EvidenceFactory::availableInfoFields, Subscribed VCF INFO field: {} is not available from the VCF file", subscribed_item);

    } else {

      active_info_map_[subscribed_item] = result->second;

    }

  }
*/
  ExecEnv::log().info("Defined VCF INFO fields: {}", vcf_info_map.size());
  struct InfoTypeCount info_count;
  for (auto const& [ident, info_record] : vcf_info_map) {

    ExecEnv::log().info("ID: {}, Description: {}, Type: {}, Number: {}", ident, info_record.description, info_record.type, info_record.number);
    info_count.incrementCount(InfoTypeCount::convertVCFType(info_record));

  }

  ExecEnv::log().info("Info type count Float {}", info_count.floatCount());
  ExecEnv::log().info("Info type count Integer {}", info_count.integerCount());
  ExecEnv::log().info("Info type count String {}", info_count.stringCount());
  ExecEnv::log().info("Info type count FloatArray {}", info_count.floatArrayCount());
  ExecEnv::log().info("Info type count IntegerArray {}", info_count.integerArrayCount());
  ExecEnv::log().info("Info type count StringArray {}", info_count.stringArrayCount());
  ExecEnv::log().info("Info type count FloatAllele {}", info_count.floatAlleleCount());
  ExecEnv::log().info("Info type count IntegerAllele {}", info_count.integerAlleleCount());
  ExecEnv::log().info("Info type count StringAllele {}", info_count.stringAlleleCount());
  ExecEnv::log().info("Info type count Boolean {}", info_count.booleanCount());


  if (info_evidence_header_->setupEvidenceHeader(active_info_map_, info_evidence_header_)) {

    ExecEnv::log().info("Subscribing to: {} VCF INFO fields", active_info_map_.size());

  } else {

    ExecEnv::log().error("EvidenceFactory::availableInfoFields, Problem Subscribing to: {} VCF INFO fields", active_info_map_.size());

  }


}
