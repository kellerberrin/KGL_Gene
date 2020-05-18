//
// Created by kellerberrin on 10/5/20.
//

#include "kgl_variant_factory_vcf_parse_info.h"

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

            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

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

            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

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

      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

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

      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

    }

  }

  return true;

}

// Slightly slower but more general info parser.
bool kgl::VCFInfoParser::infoArrayParser() {

  enum class ParserStates { KeyToken, ValueToken} parser_state = ParserStates::KeyToken;
  size_t info_length = info_view_.length();
  size_t key_token_count = 0;
  size_t key_token_offset = 0;
  size_t value_token_count = 0;
  size_t value_token_offset = 0;
  std::string_view key_view;
  std::vector<std::string_view> value_vector;
  const std::vector<std::string_view> empty_view_vector;
  const std::string_view empty_view;

  for (size_t index = 0; index < info_length; ++index) {

    switch(parser_state) {

      case ParserStates::KeyToken:
        if (info_view_[index] == INFO_VALUE_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          parser_state = ParserStates::ValueToken;
          value_token_offset = index + 1;
          value_token_count = 0;
          value_vector.clear();

        } else if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          key_view = info_view_.substr(key_token_offset, key_token_count);
          auto result = parsed_array_map_.emplace(key_view, empty_view_vector);
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <vector> pair, (duplicate)");

          }
          key_token_offset = index + 1;
          key_token_count = 0;

        } else {

          ++key_token_count;

        }
        break;

      case ParserStates::ValueToken:
        if (info_view_[index] == INFO_FIELD_DELIMITER_) {

          value_vector.emplace_back(info_view_.substr(value_token_offset, value_token_count));
          auto result = parsed_array_map_.emplace(key_view, value_vector);
          if (not result.second) {

            ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

          }
          parser_state = ParserStates::KeyToken;
          key_token_offset = index + 1;
          key_token_count = 0;

        } else if (info_view_[index] == INFO_VECTOR_DELIMITER_) {

          value_vector.emplace_back(info_view_.substr(value_token_offset, value_token_count));
          value_token_offset = index + 1;
          value_token_count = 0;

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

      ExecEnv::log().error("VCFInfoParser::infoTokenParser, Final_Parser_State=ValueToken, Final Value Token Offset: {}, Size: {} not equal the Info size : {}",
                           value_token_offset, value_token_count, info_length);

      return false;

    }

    if (value_token_count == 0) {

      value_vector.emplace_back(empty_view);

    } else {

      value_vector.emplace_back(info_view_.substr(value_token_offset, value_token_count));

    }
    auto result = parsed_array_map_.emplace(key_view, value_vector);
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

    }

  } else {

    // Check key_token_offset + key_token_count equals the length of the info field.
    // and all characters have been consumed.
    if (key_token_offset + key_token_count != info_length) {

      ExecEnv::log().error("VCFInfoParser::infoTokenParser, Final_Parser_State=KeyToken, Final Key Token Offset: {}, Size: {} not equal the Info size : {}",
                           key_token_offset, key_token_count, info_length);

      return false;

    }

    auto result = parsed_array_map_.emplace(info_view_.substr(key_token_offset, key_token_count), empty_view_vector);
    if (not result.second) {

      ExecEnv::log().warn("VCFInfoParser::infoTokenParser, cannot insert <key>, <value> pair, (duplicate)");

    }

  }

  return true;

}



// These three functions are parser specific

#ifdef USE_TOKEN_PARSER

bool kgl::VCFInfoParser::getInfoBoolean(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    if (key_it->second.second != 0) {

      ExecEnv::log().warn("VCFInfoParser::getInfoBoolean, info field key :{} expected 0 values,  size: {}", key, key_it->second.second);

    }
    return true;

  }

  return false;

}


kgl::InfoParserString kgl::VCFInfoParser::getInfoString(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    if (key_it->second.second != 1) {

      ExecEnv::log().warn("VCFInfoParser::getInfoString, info field key: {} expected 1 value, size:{}", key, key_it->second.second);

      if (key_it->second.second != 0) {

        return std::string();

      } else {

        std::vector<std::string> string_vector =Utility::char_tokenizer(std::string(key_it->second.first), INFO_VECTOR_DELIMITER_);
        return string_vector.front();

      }

    }

    return std::string(key_it->second.first);

  }

  return std::nullopt;

}



kgl::InfoParserStringArray kgl::VCFInfoParser::getInfoStringArray(const std::string& key) const {

  auto key_it = parsed_token_map_.find(key);
  if (key_it != parsed_token_map_.end()) {

    if (key_it->second.second == 1) {

      return std::vector<std::string>{ std::string(key_it->second.first) };

    }

    return Utility::char_tokenizer(std::string(key_it->second.first), INFO_VECTOR_DELIMITER_);

  }

  return std::nullopt;

}

#else


bool kgl::VCFInfoParser::getInfoBoolean(const std::string& key) const {

  auto key_it = parsed_array_map_.find(key);
  if (key_it != parsed_array_map_.end()) {

    if (not key_it->second.empty()) {

      ExecEnv::log().warn("VCFInfoParser::getInfoBoolean, info field key :{} expected 0 values,  size: {}", key, key_it->second.size());

    }
    return true;

  }

  return false;

}


kgl::InfoParserString kgl::VCFInfoParser::getInfoString(const std::string& key) const {

  auto key_it = parsed_array_map_.find(key);
  if (key_it != parsed_array_map_.end()) {

    if (key_it->second.size() != 1) {

      ExecEnv::log().warn("VCFInfoParser::getInfoString, info field key: {} expected 1 value, size:{}", key, key_it->second.size());

    }

    if (not key_it->second.empty()) {

      return std::string(key_it->second[0]);

    } else {

      return std::string();  // Return the empty string.

    }

  }

  return std::nullopt;

}



kgl::InfoParserStringArray kgl::VCFInfoParser::getInfoStringArray(const std::string& key) const {

  auto key_it = parsed_array_map_.find(key);
  if (key_it != parsed_array_map_.end()) {

    std::vector<std::string> string_vector;  //pre-allocate for efficiency
    string_vector.reserve(key_it->second.size());
    for (auto const& item : key_it->second) {

      string_vector.emplace_back(std::string(item));

    }

    return string_vector;

  }

  return std::nullopt;

}


#endif


kgl::InfoParserInteger kgl::VCFInfoParser::getInfoInteger(const std::string& key) const {

  InfoParserString string_value(getInfoString(key));

  if (not string_value) return std::nullopt;

  return convertToInteger(key, string_value.value());

}


kgl::InfoParserIntegerArray kgl::VCFInfoParser::getInfoIntegerArray(const std::string& key) const {

  InfoParserStringArray string_array(getInfoStringArray(key));

  if (not string_array) return std::nullopt;

  std::vector<std::optional<int64_t>> integer_vector{string_array.value().size()};   // Preallocate vector size.
  for (auto& value : string_array.value()) {

    integer_vector.push_back(convertToInteger(key, value));

  }

  return(std::move(integer_vector));

}


kgl::InfoParserFloat kgl::VCFInfoParser::getInfoFloat(const std::string& key) const {

  InfoParserString string_value(getInfoString(key));

  if (not string_value) return std::nullopt;

  return convertToFloat(key, string_value.value());

}


kgl::InfoParserFloatArray kgl::VCFInfoParser::getInfoFloatArray(const std::string& key) const {

  InfoParserStringArray string_array(getInfoStringArray(key));

  if (not string_array) return std::nullopt;

  std::vector<std::optional<double>> float_vector{string_array.value().size()};   // Preallocate vector size.
  for (auto& value : string_array.value()) {

    float_vector.push_back(convertToFloat(key, value));

  }

  return(std::move(float_vector));

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


bool kgl::VCFInfoParser::compareParsers() {
  static std::mutex func_lock;
  std::scoped_lock<std::mutex> raii(func_lock);

  bool bool_result = true;
  if (parsed_token_map_.size() != parsed_array_map_.size()) {

    ExecEnv::log().error("VCFInfoParser::compareParsers, Parser Size Mismatch, Token map size {}, Array map size {}", parsed_token_map_.size(), parsed_array_map_.size());
    bool_result = false;

  }

  for (auto const& [key, value] : parsed_token_map_) {

    auto result = parsed_array_map_.find(key);

    if (result == parsed_array_map_.end()) {

      ExecEnv::log().error("VCFInfoParser::compareParsers, Token map key {} not found in Array map", std::string(key));
      bool_result = false;

    } else {

      if (result->second.size() != value.second) {

        ExecEnv::log().error("VCFInfoParser::compareParsers, Token map key: {} value size {} not equal to Array map vector size: {}\nfield value {}",
                             std::string(key), value.second, result->second.size(), std::string(value.first));
        bool_result = false;

      }

      if (result->second.size() == 1) {

        if (std::string(result->second.front()) != std::string(value.first)) {

          ExecEnv::log().error("VCFInfoParser::compareParsers, Key: {}, Token map value: {} different to Array map value: {}",
                                std::string(key), std::string(value.first), std::string(result->second.front()));

          bool_result = false;

        }

      } else {

        std::vector<std::string> string_vector = Utility::char_tokenizer(std::string(value.first), INFO_VECTOR_DELIMITER_);

        std::vector<std::string> array_string_vector;
        array_string_vector.reserve(result->second.size());
        for (auto const& vector_value : result->second) {

          array_string_vector.emplace_back(std::string(vector_value));

        }

        if (string_vector.size() != array_string_vector.size()) {

          ExecEnv::log().error("VCFInfoParser::compareParsers, Key: {}, Token vector size: {} different to Array vector size: {}",
                               std::string(key), string_vector.size(), array_string_vector.size());

          for (size_t index = 0; index < string_vector.size(); ++index) {

            ExecEnv::log().error("VCFInfoParser::compareParsers, Key: {}, Token vector[{}] = {}", std::string(key), index, string_vector[index]);

          }
          for (size_t index = 0; index < array_string_vector.size(); ++index) {

            ExecEnv::log().error("VCFInfoParser::compareParsers, Key: {}, Array vector[{}] = {}", std::string(key), index, array_string_vector[index]);

          }


          bool_result = false;

        } else {

          for (size_t index = 0; index < string_vector.size(); ++index) {

            if (string_vector[index] != array_string_vector[index]) {

              ExecEnv::log().error("VCFInfoParser::compareParsers, Key: {}, Token vector[{}] = {}  different to Array vector[{}]={}",
                                   std::string(key), index, string_vector[index] , index, array_string_vector[index]);

              bool_result = false;

            }

          }

        }

      }

    }

  }

  return bool_result;

}



