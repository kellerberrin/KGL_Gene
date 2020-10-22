//
// Created by kellerberrin on 10/5/20.
//

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kel_utility.h"

#include <string_view>

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Efficient parser for the info field.
// Implemented as a simple finite state parser.
bool kgl::VCFInfoParser::infoTokenParser() {

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
          auto result = parsed_token_map_.try_emplace(key_view, InfoParserToken(empty_value_view, value_sub_field_count));
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 1. cannot insert <key> : '{}', <value> pair, (duplicate)", std::string(key_view));
            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 1. Info : {} ", std::string(info_view_));

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
          auto result = parsed_token_map_.try_emplace(key_view, InfoParserToken(value_view, value_sub_field_count));
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 2. cannot insert <key>: '{}', <value> pair, (duplicate)", std::string(key_view));
            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 2. Info : {} ", std::string(info_view_));

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

      ExecEnv::log().error("VCFInfoParser::infoTokenParser, Final_Parser_State=ValueToken, Final Value Token Offset: {}, Size: {} not equal the Info size : {}",
                           value_token_offset, value_token_count, info_length);

      return false;

    }

    value_view = info_view_.substr(value_token_offset, value_token_count);
    ++value_sub_field_count;
    auto result = parsed_token_map_.try_emplace(key_view, InfoParserToken(value_view, value_sub_field_count));
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 3. cannot insert <key> : '{}', <value> pair, (duplicate)", std::string(value_view));
      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 3. Info : {} ", std::string(info_view_));

    }

  } else {

    // Check key_token_offset + key_token_count equals the length of the info field.
    // and all characters have been consumed.
    if (key_token_offset + key_token_count != info_length) {

      ExecEnv::log().error("VCFInfoParser::infoTokenParser, Final_Parser_State=KeyToken, Final Key Token Offset: {}, Size: {} not equal the Info size : {}",
                           key_token_offset, key_token_count, info_length);

      return false;

    }

    key_view = info_view_.substr(key_token_offset, key_token_count);
    value_sub_field_count = 0;
    auto result = parsed_token_map_.try_emplace(key_view, InfoParserToken(empty_value_view, value_sub_field_count));
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 4. cannot insert <key> : '{}', <value> pair, (duplicate)", std::string(key_view));
      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, 4. Info : {} ", std::string(info_view_));

    }

  }

  return true;

}


// These three functions are parser specific


[[nodiscard]] std::optional<kgl::InfoParserToken> kgl::VCFInfoParser::getToken(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    return key_it->second;

  }

  return std::nullopt;

}


kgl::InfoIntegerType kgl::VCFInfoParser::convertToInteger(const std::string& value) {

  try {

#ifdef USE_64BIT_TYPES

    return std::stoll(std::string(value));

#else

    return std::stol(std::string(value));

#endif

  }
  catch(std::out_of_range& e) {

    ExecEnv::log().warn("VCFInfoParser::convertToInteger, Exception:Out of of Range,  Value: {}", value);
    return MISSING_VALUE_INTEGER_;

  }
  catch(std::invalid_argument& e) {

    if (std::string(value) == INFO_VECTOR_MISSING_VALUE_STR_) {

      return MISSING_VALUE_INTEGER_;

    }

    ExecEnv::log().error("VCFInfoParser::convertToInteger, Exception:Invalid Argument, Value {}", value);
    return MISSING_VALUE_INTEGER_;

  }
  catch(std::exception& e) {

    ExecEnv::log().error("VCFInfoParser::convertToInteger, Exception:Unknown, Value {}", value);
    return MISSING_VALUE_INTEGER_;

  }

}


kgl::InfoFloatType kgl::VCFInfoParser::convertToFloat(const std::string& value) {

  // Treat text 'NaN's as missing values.
  if (value.size() == 3) {

    const std::string uc_nan{"NAN"};
    if (Utility::toupper(value) == uc_nan) {

      return MISSING_VALUE_FLOAT_;

    }

  }

  try {

#ifdef USE_64BIT_TYPES

    return std::stod(value);

#else


    return std::stof(value);

#endif

  }
  catch(std::out_of_range& e) {

    const std::string underflow_exponent{"E-"};
    const std::string overflow_exponent{"E"};
    std::string uc_value = Utility::toupper(value);

    if (uc_value.find(underflow_exponent) != std::string::npos) {

      // Negative exponent exists and therefore an underflow.
      return std::numeric_limits<InfoFloatType>::min();

    } else if (uc_value.find(overflow_exponent) != std::string::npos) {

      // Exponent exists so assume an overflow.
      return std::numeric_limits<InfoFloatType>::max();

    } else { // Another floating format, rather than test exhaustively, just give up and return a non-value.

      ExecEnv::log().warn("VCFInfoParser::convertToFloat, Exception::Unknown Out of Range Error,  Value: {}", value);
      return MISSING_VALUE_FLOAT_;

    }

  }
  catch(std::invalid_argument& e) {

    if (std::string(value) == INFO_VECTOR_MISSING_VALUE_STR_) {

      return MISSING_VALUE_FLOAT_;

    }

    ExecEnv::log().error("VCFInfoParser::convertToFloat, Exception:Invalid Argument, Value {}", value);
    return MISSING_VALUE_FLOAT_;

  }
  catch(std::exception& e) {

    ExecEnv::log().error("VCFInfoParser::convertToFloat, Exception:Unknown, Value {}", value);
    return MISSING_VALUE_FLOAT_;

  }

}



