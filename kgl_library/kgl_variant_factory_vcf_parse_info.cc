//
// Created by kellerberrin on 10/5/20.
//

#include "kgl_variant_factory_vcf_parse_info.h"


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Efficient parser for the info field.
// Implemented as a simple finite state parser.
bool kgl::VCFInfoParser::parseInfo() {

  enum class ParserStates { KeyToken, ValueToken} parser_state = ParserStates::KeyToken;
  size_t info_length = info_view_.length();
  size_t key_token_count = 0;
  size_t key_token_offset = 0;
  size_t value_token_count = 0;
  size_t value_token_offset = 0;
  std::string_view key_view;
  std::string_view value_view;
  const std::string_view empty_value_view;

  for (size_t index = 0; index < info_length; ++index) {

    switch(parser_state) {

      case ParserStates::KeyToken:
        if (info_view_[index] == INFO_VALUE_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          parser_state = ParserStates::ValueToken;
          value_token_offset = index + 1;
          value_token_count = 0;

        } else if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          auto result = parsed_token_map_.emplace(key_view, empty_value_view);
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

          }
          key_token_offset = index + 1;
          key_token_count = 0;

        } else {

          ++key_token_count;

        }
        break;

      case ParserStates::ValueToken:
        if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          value_view = info_view_.substr(value_token_offset, value_token_count);
          auto result = parsed_token_map_.emplace(key_view, value_view);
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

          }
          parser_state = ParserStates::KeyToken;
          key_token_offset = index + 1;
          key_token_count = 0;

        } else {

          ++value_token_count;

        }
        break;

    } // switch

  } // for loop

  // Process final token.
  if (parser_state == ParserStates::ValueToken) {

    // Check value_token_offset + value_token_count equals the length of the info field.
    // and all characters have been consumed.
    if (value_token_offset + value_token_count != info_length) {

      ExecEnv::log().error("VCFInfoParser::parseInfo, Final_Parser_State=ValueToken, Final Value Token Offset: {}, Size: {} not equal the Info size : {}",
                           value_token_offset, value_token_count, info_length);

      return false;

    }

    value_view = info_view_.substr(value_token_offset, value_token_count);
    auto result = parsed_token_map_.emplace(key_view, value_view);
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

    }

  } else {

    // Check key_token_offset + key_token_count equals the length of the info field.
    // and all characters have been consumed.
    if (key_token_offset + key_token_count != info_length) {

      ExecEnv::log().error("VCFInfoParser::parseInfo, Final_Parser_State=KeyToken, Final Key Token Offset: {}, Size: {} not equal the Info size : {}",
                           key_token_offset, key_token_count, info_length);

      return false;

    }

    key_view = info_view_.substr(key_token_offset, key_token_count);
    auto result = parsed_token_map_.emplace(key_view, empty_value_view);
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

    }

  }

  return true;

}


std::optional<std::string> kgl::VCFInfoParser::getInfoField(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    return std::string(key_it->second);

  }

  return std::nullopt;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// An indexed map of InfoEvidenceIndex. There is only one of these held by all variants with INFO evidence fields.

bool kgl::InfoEvidenceHeader::setupEvidenceHeader(const VCFInfoRecordMap& vcf_info_map, std::shared_ptr<const InfoEvidenceHeader> self_ptr) {

  // Create a map of indexes.
  info_index_map_.clear();
  for (auto const& [ident , info_record] : vcf_info_map) {

    InfoEvidenceType type = convertVCFType(info_record);
    InfoEvidenceIndex info_index(info_record, type, 0, self_ptr);
    auto result = info_index_map_.emplace(ident, info_index);
    if (not result.second) {

      ExecEnv::log().error("InfoEvidenceHeader::setupEvidenceHeader, could not add index for info field: {}", ident);

    }

  }

  return true;
}

enum class InfoEvidenceType { Float, Integer, String, FloatArray, IntegerArray, StringArray, Boolean, NotImplemented};

kgl::InfoEvidenceType kgl::InfoEvidenceHeader::convertVCFType(const VCFInfoRecord& vcf_info_item) {

  if (vcf_info_item.type == INTEGER_) {

    if (vcf_info_item.number == SCALAR_) {

      return InfoEvidenceType::Integer;

    } else {

      return InfoEvidenceType::IntegerArray;

    }

  } else if (vcf_info_item.type == FLOAT_) {

    if (vcf_info_item.number == SCALAR_) {

      return InfoEvidenceType::Float;

    } else {

      return InfoEvidenceType::FloatArray;

    }

  } else if (vcf_info_item.type == FLAG_) {

    if (vcf_info_item.number == FLAG_SCALAR_) {

      return InfoEvidenceType::Boolean;

    } else {

      ExecEnv::log().warn("InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
                          vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
      return InfoEvidenceType::NotImplemented;

    }

  } else if (vcf_info_item.type == CHAR_STRING_ or vcf_info_item.type == STRING_) {

    if (vcf_info_item.number == SCALAR_) {

      return InfoEvidenceType::String;

    } else {

      return InfoEvidenceType::StringArray;

    }

  } else {

    ExecEnv::log().warn("InfoEvidenceHeader::convertVCFType, Ident: {}, Description: {} , Info Type: {}, Number: {} not implemented",
                         vcf_info_item.ID, vcf_info_item.description, vcf_info_item.type, vcf_info_item.number);
    return InfoEvidenceType::NotImplemented;

  }

}


std::optional<const kgl::InfoEvidenceIndex> kgl::InfoEvidenceHeader::getIndex(const std::string&) const {

  return std::nullopt;

}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The evidence factory creates a common evidence lookup object for all variants (to optimize memory usage).
// The evidence factory also creates an evidence object for each variant (data only).


std::optional<std::unique_ptr<kgl::InfoDataBlock>> kgl::EvidenceFactory::createVariantEvidence(std::string&& info) {

  // If no Info fields have been requested, then just return std::nullopt
  if (evidence_map_.empty()) {

    return std::nullopt;

  }

  // Else fill up a data block and return it.
  if (not active_info_map_.empty()) {

    // Parse the info line.
    VCFInfoParser info_parser_(std::move(info));

    for (auto const& subscribed_info : active_info_map_) {

      auto vcf_info_item = info_parser_.getInfoField(subscribed_info.first);

      if (vcf_info_item) {

        // Mark the bitfield as missing

      } else {

        // Mark the bitfield as present.

      }

    }

  }

  return std::make_unique<InfoDataBlock>(info_evidence_header_);

}


void kgl::EvidenceFactory::availableInfoFields(const VCFInfoRecordMap& vcf_info_map) {

  for (auto const& subscribed_item : evidence_map_) {

    auto result = vcf_info_map.find(subscribed_item);

    if (result == vcf_info_map.end()) {

      ExecEnv::log().warn("EvidenceFactory::availableInfoFields, Subscribed VCF INFO field: {} is not available from the VCF file", subscribed_item);

    } else {

      active_info_map_[subscribed_item] = result->second;

    }

  }

  ExecEnv::log().info("Defined VCF INFO fields: {}", vcf_info_map.size());
  for (auto const& [ident, info_record] : vcf_info_map) {

    InfoEvidenceHeader::convertVCFType(info_record);

  }


  if (info_evidence_header_->setupEvidenceHeader(active_info_map_, info_evidence_header_)) {

    ExecEnv::log().info("Subscribing to: {} VCF INFO fields", active_info_map_.size());

  } else {

    ExecEnv::log().error("EvidenceFactory::availableInfoFields, Problem Subscribing to: {} VCF INFO fields", active_info_map_.size());

  }


}

