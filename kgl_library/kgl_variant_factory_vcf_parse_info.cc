//
// Created by kellerberrin on 10/5/20.
//

#include "kgl_variant_factory_vcf_parse_info.h"

#include <string_view>

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
  size_t value_sub_field_count = 0;  // Count the number of sub fields in a value e.g. "9,9,9" = 3.
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
          value_sub_field_count = 0;

        } else if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          value_sub_field_count = 0;
          auto result = parsed_token_map_.emplace(key_view, InfoParserToken(empty_value_view, value_sub_field_count));
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
          ++value_sub_field_count;
          auto result = parsed_token_map_.emplace(key_view, InfoParserToken(value_view, value_sub_field_count));
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

          }
          parser_state = ParserStates::KeyToken;
          key_token_offset = index + 1;
          key_token_count = 0;

        } else {

          if (info_view_[index] == INFO_VECTOR_DELIMITER_) ++value_sub_field_count;
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
    ++value_sub_field_count;
    auto result = parsed_token_map_.emplace(key_view, InfoParserToken(value_view, value_sub_field_count));
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
    value_sub_field_count = 0;
    auto result = parsed_token_map_.emplace(key_view, InfoParserToken(empty_value_view, value_sub_field_count));
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::parseInfo, cannot insert <key>, <value> pair, (duplicate)");

    }

  }

  return true;

}


std::optional<std::string> kgl::VCFInfoParser::getInfoString(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    return std::string(key_it->second.first);

  }

  return std::nullopt;

}


std::optional<std::vector<std::string>> kgl::VCFInfoParser::getInfoStringArray(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    return Utility::tokenizer(std::string(key_it->second.first), INFO_VECTOR_DELIMITER_STR_);

  }

  return std::nullopt;

}

std::optional<int64_t> kgl::VCFInfoParser::getInfoInteger(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    return convertToInteger(key, std::string(key_it->second.first));

  }

  return std::nullopt; // field not found.

}


std::optional<std::vector<std::optional<int64_t>>> kgl::VCFInfoParser::getInfoIntegerArray(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    std::vector<std::string> value_array = Utility::tokenizer(std::string(key_it->second.first), INFO_VECTOR_DELIMITER_STR_);
    std::vector<std::optional<int64_t>> integer_vector{value_array.size()};   // Preallocate vector size.
    for (auto& value : value_array) {

      integer_vector.push_back(convertToInteger(key, value));

    }

    return(std::move(integer_vector));

  }

  return std::nullopt; // field not found.

}


std::optional<double> kgl::VCFInfoParser::getInfoFloat(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    return convertToFloat(key, std::string(key_it->second.first));

  }

  return std::nullopt; // field not found.

}


std::optional<std::vector<std::optional<double>>> kgl::VCFInfoParser::getInfoFloatArray(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    std::vector<std::string> value_array = Utility::tokenizer(std::string(key_it->second.first), INFO_VECTOR_DELIMITER_STR_);
    std::vector<std::optional<double>> float_vector{value_array.size()};   // Preallocate vector size.
    for (auto& value : value_array) {

      float_vector.push_back(convertToFloat(key, value));

    }

    return(std::move(float_vector));

  }

  return std::nullopt; // field not found.

}


std::optional<int64_t> kgl::VCFInfoParser::convertToInteger(const std::string& key, const std::string& value) {

  try {

    return std::stoll(std::string(value));

  }
  catch(std::out_of_range& e) {

    ExecEnv::log().warn("VCFInfoParser::convertToInteger, Exception:Out of of Range,  Field:{}, Value: {}", key, value);
    return std::nullopt;

  }
  catch(std::invalid_argument& e) {

    if (std::string(value) == INFO_VECTOR_MISSING_VALUE_STR_) {

      return std::nullopt;

    }

    ExecEnv::log().error("VCFInfoParser::convertToInteger, Exception:Invalid Argument, Field {}, Value {}",key, value);
    return std::nullopt;

  }
  catch(std::exception& e) {

    ExecEnv::log().error("VCFInfoParser::convertToInteger, Exception:Unknown, Field {}, Value {}",key, value);
    return std::nullopt;

  }

}


std::optional<double> kgl::VCFInfoParser::convertToFloat(const std::string& key, const std::string& value) {

  try {

    return std::stod(std::string(value));

  }
  catch(std::out_of_range& e) {

    ExecEnv::log().warn("VCFInfoParser::convertToFloat, Exception:Out of of Range,  Field:{}, Value: {}", key, value);
    return std::nullopt;

  }
  catch(std::invalid_argument& e) {

    if (std::string(value) == INFO_VECTOR_MISSING_VALUE_STR_) {

      return std::nullopt;

    }

    ExecEnv::log().error("VCFInfoParser::convertToFloat, Exception:Invalid Argument, Field {}, Value {}",key, value);
    return std::nullopt;

  }
  catch(std::exception& e) {

    ExecEnv::log().error("VCFInfoParser::convertToFloat, Exception:Unknown, Field {}, Value {}",key, value);
    return std::nullopt;

  }

}


