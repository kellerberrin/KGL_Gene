//
// Created by kellerberrin on 23/5/20.
//

#include "kgl_variant_factory_vcf_evidence_memory.h"
#include "kgl_variant_factory_vcf_evidence.h"



namespace kgl = kellerberrin::genome;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// data_item_count is the number of data items, for a float array of size 3 this is 1 for 1 array. For a string array of size 3 it is 3.
// data_item_size is the number of underlying data items. For float array of size 3 this is 3.
// For a string array this will be the total number of characters in all strings. For a string array of size 3 data_item_count=3, but data_item_size is size of
// all three strings added together. The total characters required for all 3 strings.
kgl::ItemOffset kgl::DataInfoTypeCount::staticIncrementAndAllocate(InfoEvidenceIntern internal_type) {

  ItemOffset item_offset;
  switch (internal_type) {

    case InfoEvidenceIntern::intern_char:  {

      item_offset.offset = char_count_;
      item_offset.is_array = false;
      ++char_count_;

    }// Size is fixed (boolean variables) and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_unity_integer_array:
    case InfoEvidenceIntern::intern_integer: {

      item_offset.offset = integer_count_;
      item_offset.is_array = false;
      ++integer_count_;

    }// Size is fixed and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_unity_float_array:
    case InfoEvidenceIntern::intern_float: {

      item_offset.offset = float_count_;
      item_offset.is_array = false;
      ++float_count_;  // should be 1

    }// Size is fixed and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_string: {

      item_offset.offset = string_count_;
      item_offset.is_array = false;
      ++string_count_;   // allocate a std::string_view

    } // Size varies between records.
      return item_offset;

    case InfoEvidenceIntern::intern_integer_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    } // Size is fixed and known at Info data subscription time.
      return item_offset;

    case InfoEvidenceIntern::intern_float_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    }
      return item_offset;


    case InfoEvidenceIntern::intern_string_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    } // Size varies between records.
      return item_offset;

    case InfoEvidenceIntern::intern_unity_string_array: {

      item_offset.offset = array_count_;
      item_offset.is_array = true;
      ++array_count_;

    }   // Size varies between records.
      return item_offset;

    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
    default:
      return item_offset;

  }

}


bool kgl::DataInfoTypeCount::dynamicIncrementAndAllocate(const InfoSubscribedField& subscribed_field, const InfoParserToken& token) {

  InfoEvidenceIntern internal_type = subscribed_field.evidenceType().InternalInfoType();

  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      if (token.second != 0) {

        ExecEnv::log().warn("DataInfoTypeCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
        , std::string(token.first), token.second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
        return false;
      }
      break;

    case InfoEvidenceIntern::intern_integer:
    case InfoEvidenceIntern::intern_float:
      if (token.second != 1) {

        ExecEnv::log().warn("DataInfoTypeCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
        , std::string(token.first), token.second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
        return false;
      }
      break;

    case InfoEvidenceIntern::intern_unity_integer_array:
      if (token.second > 1) {

        ++unity_array_count_;
        integer_count_ += token.second;

      }
      break;

    case InfoEvidenceIntern::intern_unity_float_array:
      if (token.second > 1) {

        ++unity_array_count_;
        float_count_ += token.second;

      }
      break;

    case InfoEvidenceIntern::intern_string: {

      char_count_ += token.first.size(); // total char size of all strings.

    }
      break;

    case InfoEvidenceIntern::intern_integer_array: {

      integer_count_ += token.second; // size of the array

    } // Size is fixed and known at Info data subscription time.
      break;

    case InfoEvidenceIntern::intern_float_array: {

      float_count_ += token.second; // size of the array

    } // Size is fixed and known at Info data subscription time.
      break;

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array: {

      char_count_ += token.first.size() - (token.second - 1); // total char size of all strings, less the delimiter chars.
      string_count_ += token.second;   // number strings, allocate a vector of std::string_views

    } // Size varies between records.
      break;

    case InfoEvidenceIntern::NotImplemented:  // Trivially fixed.
    default:
      break;

  }

  return true;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::InfoDataBlock::indexAndVerify( const InfoSubscribedField& subscribed_field,     // The field index
                                         const std::optional<InfoParserToken>& token, // The parsed token to be placed in memory
                                         DataInfoTypeCount& dynamic_accounting,  // accounting objects to check index, data and memory integrity
                                         DataInfoTypeCount& static_accounting,
                                         std::vector<ItemOffset>& index_accounting) {

  InfoEvidenceIntern internal_type = subscribed_field.evidenceType().InternalInfoType();

  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      return internChar(subscribed_field, token, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_integer:
      return internInteger(subscribed_field, token, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_float:
      return internFloat(subscribed_field, token, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_unity_integer_array:
      return internUnityInteger(subscribed_field, token, dynamic_accounting, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_unity_float_array:
      return internUnityFloat(subscribed_field, token, dynamic_accounting, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_string:
      return internString(subscribed_field, token, dynamic_accounting, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_integer_array:
      return internIntegerArray(subscribed_field, token, dynamic_accounting, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_float_array:
      return internFloatArray(subscribed_field, token, dynamic_accounting, static_accounting, index_accounting);

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array:
      return internStringArray(subscribed_field, token, dynamic_accounting, static_accounting, index_accounting);

    case InfoEvidenceIntern::NotImplemented:  // unknown internal type.
    default:
      ExecEnv::log().warn( "InfoDataBlock::indexAndVerify, Internal data type unknown for field ID:{}, Number:{}, Type:{}",
                           subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
      return false;

  }

}


// These are just function blocks that break up the unsightly mass of code in function InfoDataBlock::indexAndVerify.
bool kgl::InfoDataBlock::internChar( const InfoSubscribedField& subscribed_field,     // The field index
                                     const std::optional<InfoParserToken>& token, // The parsed token to be placed in memory
                                     DataInfoTypeCount& static_accounting,
                                     std::vector<ItemOffset>& index_accounting) {

  // check field index
  if (subscribed_field.dataOffset().offset >= type_count_.charCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} exceeds allocated data size: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, type_count_.charCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Add the data.
  char_memory_[subscribed_field.dataOffset().offset] = token ? true : false;

  return true;

}

bool kgl::InfoDataBlock::internInteger( const InfoSubscribedField& subscribed_field,
                                        const std::optional<InfoParserToken>& token,
                                        DataInfoTypeCount& static_accounting,
                                        std::vector<ItemOffset>& index_accounting) {

  if (subscribed_field.dataOffset().offset >= type_count_.integerCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} exceeds allocated data size: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, type_count_.integerCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().critical(
      "InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}, field ID:{}, Number:{}, Type:{}",
      token.value().second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
      subscribed_field.infoRecord().type);

    }

    integer_memory_[subscribed_field.dataOffset().offset] = VCFInfoParser::convertToInteger(
    std::string(token.value().first));

  } else {

    integer_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_INTEGER;

  }

  return true;

}

bool kgl::InfoDataBlock::internFloat( const InfoSubscribedField& subscribed_field,
                                      const std::optional<InfoParserToken>& token,
                                      DataInfoTypeCount& static_accounting,
                                      std::vector<ItemOffset>& index_accounting) {

  if (subscribed_field.dataOffset().offset >= type_count_.floatCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} exceeds allocated data size: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, type_count_.floatCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().critical(
      "InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}, field ID:{}, Number:{}, Type:{}",
      token.value().second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
      subscribed_field.infoRecord().type);

    }

    float_memory_[subscribed_field.dataOffset().offset] = VCFInfoParser::convertToFloat(
    std::string(token.value().first));

  } else {

    float_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_FLOAT;

  }

  return true;

}

bool kgl::InfoDataBlock::internString( const InfoSubscribedField& subscribed_field,
                                       const std::optional<InfoParserToken>& token,
                                       DataInfoTypeCount& dynamic_accounting,
                                       DataInfoTypeCount& static_accounting,
                                       std::vector<ItemOffset>& index_accounting) {

  if (subscribed_field.dataOffset().offset >= type_count_.stringCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} exceeds allocated data size: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, type_count_.floatCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().warn(
      "InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}, field ID:{}, Number:{}, Type:{}, Value: {}",
      token.value().second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
      subscribed_field.infoRecord().type, std::string(token.value().first));

    }

    char* string_start = &char_memory_[dynamic_accounting.charCount()];
    size_t string_size = token.value().first.size();
    string_memory_[subscribed_field.dataOffset().offset] = std::string_view(string_start, string_size);
    if (dynamic_accounting.charCount() + string_size > type_count_.charCount()) {

      ExecEnv::log().critical(
      "InfoDataBlock::indexAndVerify, String size+offset :{} has exceeded char vector size: {}, field ID:{}, Number:{}, Type:{}",
      dynamic_accounting.charCount() + string_size, type_count_.charCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
      subscribed_field.infoRecord().type);

    }
    std::memcpy(string_start, token.value().first.data(), string_size);
    dynamic_accounting.charCount(dynamic_accounting.charCount() + string_size);

  } else {

    string_memory_[subscribed_field.dataOffset().offset] = std::string_view();

  }

  return true;

}

bool kgl::InfoDataBlock::internIntegerArray( const InfoSubscribedField& subscribed_field,
                                             const std::optional<InfoParserToken>& token,
                                             DataInfoTypeCount& dynamic_accounting,
                                             DataInfoTypeCount& static_accounting,
                                             std::vector<ItemOffset>& index_accounting) {

  return true;

}

bool kgl::InfoDataBlock::internFloatArray( const InfoSubscribedField& subscribed_field,
                                           const std::optional<InfoParserToken>& token,
                                           DataInfoTypeCount& dynamic_accounting,
                                           DataInfoTypeCount& static_accounting,
                                           std::vector<ItemOffset>& index_accounting) {

  return true;

}

bool kgl::InfoDataBlock::internStringArray( const InfoSubscribedField& subscribed_field,
                                            const std::optional<InfoParserToken>& token,
                                            DataInfoTypeCount& dynamic_accounting,
                                            DataInfoTypeCount& static_accounting,
                                            std::vector<ItemOffset>& index_accounting) {

  return true;

}

bool kgl::InfoDataBlock::internUnityInteger( const InfoSubscribedField& subscribed_field,
                                             const std::optional<InfoParserToken>& token,
                                             DataInfoTypeCount& dynamic_accounting,
                                             DataInfoTypeCount& static_accounting,
                                             std::vector<ItemOffset>& index_accounting) {

  return true;

}

bool kgl::InfoDataBlock::internUnityFloat( const InfoSubscribedField& subscribed_field,
                                           const std::optional<InfoParserToken>& token,
                                           DataInfoTypeCount& dynamic_accounting,
                                           DataInfoTypeCount& static_accounting,
                                           std::vector<ItemOffset>& index_accounting) {

  return true;

}

