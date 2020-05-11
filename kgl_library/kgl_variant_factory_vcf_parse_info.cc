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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

void kgl::GnomadInfo_3_0::loadActive() {

  for (auto const& [field, field_info] : gnomad_3_0_map_) {

    if (field_info.process) {

      auto result = active_map_.emplace(field, field_info);

      if (not result.second) {

        ExecEnv::log().error("kgl::GnomadInfo_3_0::loadActive, Could not active INFO field : {}", field);

      }

    }

  }

}
