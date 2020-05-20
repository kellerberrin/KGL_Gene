//
// Created by kellerberrin on 16/5/20.
//

#include "kgl_variant_factory_vcf_evidence.h"
#include "kgl_variant_factory_vcf_parse_header.h"
#include "kgl_variant_factory_vcf_parse_info.h"

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoEvidenceIndex. There is only one of these held by all variants with INFO evidence fields.

bool kgl::InfoEvidenceHeader::setupEvidenceHeader(const InfoSubscribedMap& vcf_info_map, std::shared_ptr<const InfoEvidenceHeader> self_ptr) {

  // Create a map of indexes.
  info_index_map_.clear();
  for (auto const& [ident , info_record] : vcf_info_map) {

    InfoEvidenceTypePair type = InfoTypeCount::evidenceType(info_record.infoRecord());
    InfoEvidenceIndex info_index(info_record.infoRecord(), type, 0, self_ptr);
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


// Returns true if a number > 1. Used in the type definition lambdas.
bool kgl::InfoTypeCount::isVectorType(const std::string& type) {

  try {

    return std::stol(type) > 1;

  }
  catch(...) {

    return false;

  }

}


const kgl::InfoEvidenceTypePair kgl::InfoTypeCount::evidenceType(const VCFInfoRecord &vcf_info_item) {

  for (auto const& evidence_type : type_definitions_) {

    if (evidence_type.second(vcf_info_item.type, vcf_info_item.number)) {

      return evidence_type.first;

    }

  }

  ExecEnv::log().warn("InfoTypeCount::evidenceType, Info ID: {} Unable to find data type combinbation for Number: {}, Type: {}",
                       vcf_info_item.ID, vcf_info_item.number, vcf_info_item.type);

  return {InfoEvidenceApp::NotImplemented, InfoEvidenceImpl::NotImplemented };

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// InfoDataBlock

kgl::DataInfoTypeCount kgl::InfoDataBlockNaive::dataPayload() {


//  static std::mutex func_lock;
//  std::scoped_lock<std::mutex> raii(func_lock);


  size_t float_count = 0;
  size_t integer_count = 0;
  size_t char_count = 0;
  size_t array_count = 0;
  size_t string_count = 0;

  char_count = bool_data_.size();

  integer_count = integer_data_.size();

  float_count = float_data_.size();

  size_t char_bytes = 0;
  for (auto const& string_field :  string_data_) {

    ++string_count;
    char_bytes += string_field.size();

  }

  for (auto const& float_array :  float_array_data_) {

    if (float_array.size() <= 1) {

      ++float_count;

    } else {

      float_count += float_array.size();
      ++array_count;

    }

  }

  for (auto const&integer_array :  integer_array_data_) {

    if (integer_array.size() <= 1) {

      ++float_count;

    } else {

      integer_count += integer_array.size();
      ++array_count;

    }

  }

  for (auto const& string_array :  string_array_data_) {

    if (string_array_data_.size() == 1) {

      ++string_count;
      char_bytes += string_array_data_.front().size();

    } else if (string_array_data_.empty() ) {

      ++string_count;

    } else {

      ++array_count;

      for (auto const& string_item : string_array) {

        ++string_count;
        char_bytes += string_item.size();

      }

    }

  }

  DataInfoTypeCount info_count;

  info_count.integer_count_ = integer_count;
  info_count.float_count_ = float_count;
  info_count.char_count_ = char_count;
  info_count.array_count_ = array_count;
  info_count.string_count_ = string_count;

  return info_count;

}

void kgl::InfoDataBlockNaive::clearAll() {

// The stored data.
  bool_data_.clear();
  float_data_.clear();
  integer_data_.clear();
  string_data_.clear();
  float_array_data_.clear();
  integer_array_data_.clear();
  string_array_data_.clear();

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant (data only).


kgl::InfoDataEvidence kgl::EvidenceFactory::createVariantEvidence(std::string&& info) {

  // If no Info fields have been requested, then just return std::nullopt
  if (evidence_map_.empty()) {

    return std::nullopt;

  }

  std::unique_ptr<InfoDataBlockNaive> temp = parseSubscribed(std::move(info));
//  temp->shrinkToFit();
  DataInfoTypeCount type_count = temp->dataPayload();

//  return std::move(temp);

  std::unique_ptr<InfoDataBlock> data_block_ptr(std::make_unique<InfoDataBlock>(info_evidence_header_));

  data_block_ptr->allocateMemory(type_count);

  return std::move(data_block_ptr);

}


std::unique_ptr<kgl::InfoDataBlockNaive> kgl::EvidenceFactory::parseSubscribed(std::string&& info) {

  // Parse the VCF info line.
  VCFInfoParser info_parser(std::move(info));

  std::unique_ptr info_data_ptr(std::make_unique<InfoDataBlockNaive>(info_evidence_header_));
  // Fill all the subscribed Info field values.
  for (auto const& subscribed_info : active_info_map_) {

    switch(subscribed_info.second.infoType().InfoType()) {

      case InfoEvidenceImpl::Float: {
        info_data_ptr->float_data_.emplace_back(info_parser.getInfoFloat(subscribed_info.first));
      }
      break;

      case InfoEvidenceImpl::Integer: {

        info_data_ptr->integer_data_.emplace_back(info_parser.getInfoInteger(subscribed_info.first));

      }
      break;

      case InfoEvidenceImpl::String: {

        info_data_ptr->string_data_.emplace_back(info_parser.getInfoString(subscribed_info.first));

      }
      break;

      case InfoEvidenceImpl::Boolean: {

        info_data_ptr->bool_data_.emplace_back(info_parser.getInfoBoolean(subscribed_info.first));

      }
      break;

      case InfoEvidenceImpl::IntegerArray: {

        info_data_ptr->integer_array_data_.emplace_back(info_parser.getInfoIntegerArray(subscribed_info.first));

      }
      break;

      case InfoEvidenceImpl::FloatArray: {

        info_data_ptr->float_array_data_.emplace_back(info_parser.getInfoFloatArray(subscribed_info.first));

      }
      break;

      case InfoEvidenceImpl::StringArray: {

        info_data_ptr->string_array_data_.emplace_back(info_parser.getInfoStringArray(subscribed_info.first));

      }
      break;

      default:
        ExecEnv::log().error("EvidenceFactory::parseSubscribed, Unexpected INFO data type for field: {}", subscribed_info.first);
        break;

    }

  }

  return info_data_ptr;

}


void kgl::EvidenceFactory::availableInfoFields(const VCFInfoRecordMap& vcf_info_map) {


  all_available_map_ = vcf_info_map;

  for (auto const& [ident, vcf_info_field] : vcf_info_map) {

    InfoSubscribed subscribed(vcf_info_field);
    auto result = active_info_map_.try_emplace(ident, subscribed);

    if (not result.second) {

      ExecEnv::log().error("EvidenceFactory::availableInfoFields. (Duplicate) unable to subscribe to duplicate info field: {}", ident);

    }

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
  for (auto const& [ident, info_record] : vcf_info_map) {

    ExecEnv::log().info("ID: {}, Description: {}, Type: {}, Number: {}", ident, info_record.description, info_record.type, info_record.number);

  }


  if (info_evidence_header_->setupEvidenceHeader(active_info_map_, info_evidence_header_)) {

    ExecEnv::log().info("Subscribing to: {} VCF INFO fields", active_info_map_.size());

  } else {

    ExecEnv::log().error("EvidenceFactory::availableInfoFields, Problem Subscribing to: {} VCF INFO fields", active_info_map_.size());

  }


}
