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
kgl::ItemOffset kgl::InfoDataUsageCount::staticIncrementAndAllocate(InfoEvidenceIntern internal_type) {

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


bool kgl::InfoDataUsageCount::dynamicIncrementAndAllocate(const InfoSubscribedField& subscribed_field, const InfoParserToken& token) {

  InfoEvidenceIntern internal_type = subscribed_field.evidenceType().InternalInfoType();

  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      if (token.second != 0) {

        ExecEnv::log().warn("InfoDataUsageCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
        , std::string(token.first), token.second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
        return false;
      }
      break;

    case InfoEvidenceIntern::intern_integer:
    case InfoEvidenceIntern::intern_float:
      if (token.second != 1) {

        ExecEnv::log().warn("InfoDataUsageCount::dynamicIncrementAndAllocate, Bad size (expected 1) Token: {} size: {}, field ID:{}, Number:{}, Type:{}"
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


bool kgl::InfoDataUsageCount::operator==(const InfoDataUsageCount& cmp) const {

  bool result = unityArrayCount() == cmp.unityArrayCount() and
                arrayCount() == cmp.arrayCount() and
                floatCount() == cmp.floatCount() and
                integerCount() == cmp.integerCount() and
                stringCount() == cmp.stringCount() and
                charCount() == charCount();

  return result;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Actually allocates the raw memory for each Info field.
void kgl::InfoDataBlock::allocateMemory(const InfoDataUsageCount &type_count) {

  type_count_ = type_count;
  char_memory_ = std::make_unique<char[]>(type_count_.charCount());
  integer_memory_ = std::make_unique<InfoIntegerType[]>(type_count_.integerCount());
  float_memory_ = std::make_unique<InfoFloatType[]>(type_count_.floatCount());
  array_memory_ = std::make_unique<InfoArrayIndex[]>(type_count_.arrayCount());
  unity_array_memory_ = std::make_unique<InfoArrayIndex[]>(type_count_.unityArrayCount());
  string_memory_ = std::make_unique<std::string_view[]>(type_count_.stringCount());

}



bool kgl::InfoDataBlock::indexAndVerify(const InfoSubscribedField& subscribed_field,     // The field index
                                         const std::optional<InfoParserToken>& token, // The parsed token to be placed in memory
                                         InfoDataUsageCount& dynamic_accounting,  // accounting objects to check index, data and memory integrity
                                         InfoDataUsageCount& static_accounting) {

  InfoEvidenceIntern internal_type = subscribed_field.evidenceType().InternalInfoType();

  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      return internChar(subscribed_field, token, static_accounting);

    case InfoEvidenceIntern::intern_integer:
      return internInteger(subscribed_field, token, static_accounting);

    case InfoEvidenceIntern::intern_float:
      return internFloat(subscribed_field, token, static_accounting);

    case InfoEvidenceIntern::intern_unity_integer_array:
      return internUnityInteger(subscribed_field, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_unity_float_array:
      return internUnityFloat(subscribed_field, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_string:
      return internString(subscribed_field, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_integer_array:
      return internIntegerArray(subscribed_field, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_float_array:
      return internFloatArray(subscribed_field, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array:
      return internStringArray(subscribed_field, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::NotImplemented:  // unknown internal type.
    default:
      ExecEnv::log().warn( "InfoDataBlock::indexAndVerify, Internal data type unknown for field ID:{}, Number:{}, Type:{}",
                           subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);
      return false;

  }

}


// These are just function blocks that break up the unsightly mass of code in function InfoDataBlock::indexAndVerify.
bool kgl::InfoDataBlock::internChar(const InfoSubscribedField& subscribed_field,     // The field index
                                     const std::optional<InfoParserToken>& token, // The parsed token to be placed in memory
                                     InfoDataUsageCount& static_accounting) {

  // check field index
  if (subscribed_field.dataOffset().offset != static_accounting.charCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, static_accounting.charCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Increment accounting.
  static_accounting.charCountAdd(1);
  // Add the data.
  char_memory_[subscribed_field.dataOffset().offset] = token ? true : false;

  return true;

}

bool kgl::InfoDataBlock::internInteger(const InfoSubscribedField& subscribed_field,
                                       const std::optional<InfoParserToken>& token,
                                       InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.integerCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, static_accounting.integerCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }

  // Increment the accounting.
  static_accounting.integerCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}, field ID:{}, Number:{}, Type:{}",
      token.value().second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

    }

    integer_memory_[subscribed_field.dataOffset().offset] = VCFInfoParser::convertToInteger(std::string(token.value().first));

  } else {

    integer_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_INTEGER;

  }

  return true;

}

bool kgl::InfoDataBlock::internFloat(const InfoSubscribedField& subscribed_field,
                                     const std::optional<InfoParserToken>& token,
                                     InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.floatCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, static_accounting.floatCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Increment the accounting.
  static_accounting.floatCountAdd(1);
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

bool kgl::InfoDataBlock::internString(const InfoSubscribedField& subscribed_field,
                                      const std::optional<InfoParserToken>& token,
                                      InfoDataUsageCount& dynamic_accounting,
                                      InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.stringCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, static_accounting.stringCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Update the static accounting
  static_accounting.stringCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().warn("InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}, field ID:{}, Number:{}, Type:{}, Value: {}",
      token.value().second, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type, std::string(token.value().first));

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
    // Update the dynamic accounting.
    dynamic_accounting.charCountAdd(string_size);

  } else {
// Empty string.
    string_memory_[subscribed_field.dataOffset().offset] = std::string_view();

  }

  return true;

}


bool kgl::InfoDataBlock::internIntegerArray(const InfoSubscribedField& subscribed_field,
                                            const std::optional<InfoParserToken>& token,
                                            InfoDataUsageCount& dynamic_accounting,
                                            InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.arrayCount()) {

    ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
                            subscribed_field.dataOffset().offset, static_accounting.arrayCount(), subscribed_field.infoRecord().ID,
                            subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Update the static accounting
  static_accounting.arrayCountAdd(1);
  // add the data.
  if (token) {

    size_t array_size = token.value().second;
    array_memory_[subscribed_field.dataOffset().offset].infoOffset(dynamic_accounting.integerCount());
    array_memory_[subscribed_field.dataOffset().offset].infoSize(array_size);
    array_memory_[subscribed_field.dataOffset().offset].infoVariableIndex(subscribed_field.fieldIndex());

    if (dynamic_accounting.integerCount() + array_size > type_count_.integerCount()) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Integer Array size+offset :{} has exceeded Integer vector size: {}, field ID:{}, Number:{}, Type:{}",
                              dynamic_accounting.integerCount() + array_size, type_count_.integerCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
                              subscribed_field.infoRecord().type);

    }
    // Get a vector of integers
    InfoParserIntegerArray integer_array =  VCFInfoParser::getInfoIntegerArray(std::string(token.value().first));
    if (integer_array.size() != array_size) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Integer Array size : {} differs from token size: {}, field ID:{}, Number:{}, Type:{}",
                              integer_array.size(), array_size, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

    }

    // Add to the memory array.
    size_t memory_index = dynamic_accounting.integerCount();
    for (InfoIntegerType integer_value : integer_array) {

      integer_memory_[memory_index] = integer_value;
      ++memory_index;

    }

    // Update the dynamic accounting.
    dynamic_accounting.integerCountAdd(array_size);

  } else {
// Empty vector
    array_memory_[subscribed_field.dataOffset().offset].infoOffset(0);
    array_memory_[subscribed_field.dataOffset().offset].infoSize(0);
    array_memory_[subscribed_field.dataOffset().offset].infoVariableIndex(subscribed_field.fieldIndex());

  }

  return true;

}


bool kgl::InfoDataBlock::internFloatArray(const InfoSubscribedField& subscribed_field,
                                          const std::optional<InfoParserToken>& token,
                                          InfoDataUsageCount& dynamic_accounting,
                                          InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.arrayCount()) {

    ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
                            subscribed_field.dataOffset().offset, static_accounting.arrayCount(), subscribed_field.infoRecord().ID,
                            subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Update the static accounting
  static_accounting.arrayCountAdd(1);
  // add the data.
  if (token) {

    size_t array_size = token.value().second;
    array_memory_[subscribed_field.dataOffset().offset].infoOffset(dynamic_accounting.floatCount());
    array_memory_[subscribed_field.dataOffset().offset].infoSize(array_size);
    array_memory_[subscribed_field.dataOffset().offset].infoVariableIndex(subscribed_field.fieldIndex());

    if (dynamic_accounting.floatCount() + array_size > type_count_.floatCount()) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Float Array size+offset :{} has exceeded Float vector size: {}, field ID:{}, Number:{}, Type:{}",
                              dynamic_accounting.floatCount() + array_size, type_count_.floatCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
                              subscribed_field.infoRecord().type);

    }
    // Get a vector of integers
    InfoParserFloatArray float_array =  VCFInfoParser::getInfoFloatArray(std::string(token.value().first));
    if (float_array.size() != array_size) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Float Array size : {} differs from token size: {}, field ID:{}, Number:{}, Type:{}",
                              float_array.size(), array_size, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

    }

    // Add to the memory array.
    size_t memory_index = dynamic_accounting.floatCount();
    for (InfoIntegerType float_value : float_array) {

      float_memory_[memory_index] = float_value;
      ++memory_index;

    }

    // Update the dynamic accounting.
    dynamic_accounting.floatCountAdd(array_size);

  } else {
// Empty vector
    array_memory_[subscribed_field.dataOffset().offset].infoOffset(0);
    array_memory_[subscribed_field.dataOffset().offset].infoSize(0);
    array_memory_[subscribed_field.dataOffset().offset].infoVariableIndex(subscribed_field.fieldIndex());

  }

  return true;

}


bool kgl::InfoDataBlock::internStringArray(const InfoSubscribedField& subscribed_field,
                                           const std::optional<InfoParserToken>& token,
                                           InfoDataUsageCount& dynamic_accounting,
                                           InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.arrayCount()) {

    ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
                            subscribed_field.dataOffset().offset, static_accounting.arrayCount(), subscribed_field.infoRecord().ID,
                            subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Update the static accounting
  static_accounting.arrayCountAdd(1);
  // add the data.
  if (token) {

    size_t array_size = token.value().second;
    array_memory_[subscribed_field.dataOffset().offset].infoOffset(dynamic_accounting.stringCount());
    array_memory_[subscribed_field.dataOffset().offset].infoSize(array_size);
    array_memory_[subscribed_field.dataOffset().offset].infoVariableIndex(subscribed_field.fieldIndex());

    if (dynamic_accounting.stringCount() + array_size > type_count_.stringCount()) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, String Array size+offset :{} has exceeded String vector size: {}, field ID:{}, Number:{}, Type:{}",
                              dynamic_accounting.stringCount() + array_size, type_count_.stringCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
                              subscribed_field.infoRecord().type);

    }
    // Get a vector of string_views.
    std::vector<std::string_view> string_view_array = VCFInfoParser::getInfoStringArray(token.value().first);
    if (string_view_array.size() != array_size) {

      ExecEnv::log().critical("InfoDataBlock::indexAndVerify, String Array size : {} differs from token size: {}, field ID:{}, Number:{}, Type:{}",
                              string_view_array.size(), array_size, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

    }

    // Add to the memory array.
    size_t memory_index = dynamic_accounting.stringCount();
    for (auto str_view : string_view_array) {

      char* string_address = &char_memory_[dynamic_accounting.charCount()];
      string_memory_[memory_index] = std::string_view(string_address, str_view.size());
      // update char count.
      dynamic_accounting.charCountAdd(str_view.size());

      if (dynamic_accounting.charCount() > type_count_.charCount()) {

        ExecEnv::log().critical(
        "InfoDataBlock::indexAndVerify, String size+offset :{} has exceeded char vector size: {}, field ID:{}, Number:{}, Type:{}",
        dynamic_accounting.charCount(), type_count_.charCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,subscribed_field.infoRecord().type);

      }

      std::memcpy(string_address, str_view.data(), str_view.size());
      ++memory_index;

    }

    // Update the dynamic accounting for strings.
    dynamic_accounting.stringCountAdd(array_size);

  } else {
// Empty vector
    array_memory_[subscribed_field.dataOffset().offset].infoOffset(0);
    array_memory_[subscribed_field.dataOffset().offset].infoSize(0);
    array_memory_[subscribed_field.dataOffset().offset].infoVariableIndex(subscribed_field.fieldIndex());

  }

  return true;

}

bool kgl::InfoDataBlock::internUnityInteger(const InfoSubscribedField& subscribed_field,
                                            const std::optional<InfoParserToken>& token,
                                            InfoDataUsageCount& dynamic_accounting,
                                            InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.integerCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, static_accounting.integerCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }

  // Increment the static accounting.
  static_accounting.integerCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second > 1) {
// We have a vector, which can (infrequently) occur with this internal datatype.

      // Put a missing value in the static storage allocated.
      integer_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_INTEGER;

      if (dynamic_accounting.unityArrayCount() > type_count_.unityArrayCount()) {

        ExecEnv::log().error("InfoDataBlock::indexAndVerify, Dynamic Unity Array Size: {} exceeds allocated Unity Array Size: {}",
                             dynamic_accounting.unityArrayCount(), type_count_.unityArrayCount());
        return false;
      }
      // Insert the array details into the unity array object.
      // Note the field index is not used.
      size_t array_size = token.value().second;
      unity_array_memory_[dynamic_accounting.unityArrayCount()].infoOffset(dynamic_accounting.floatCount());
      unity_array_memory_[dynamic_accounting.unityArrayCount()].infoSize(array_size);
      unity_array_memory_[dynamic_accounting.unityArrayCount()].infoVariableIndex(subscribed_field.fieldIndex());

      // Update the unity array count.
      dynamic_accounting.unityArrayCountAdd(1);

      // Copy in the integer data.
      if (dynamic_accounting.integerCount() + array_size > type_count_.integerCount()) {

        ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Integer Array size+offset :{} has exceeded Integer vector size: {}, field ID:{}, Number:{}, Type:{}",
                                dynamic_accounting.integerCount() + array_size, type_count_.integerCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
                                subscribed_field.infoRecord().type);

      }
      // Get a vector of integers
      InfoParserIntegerArray integer_array =  VCFInfoParser::getInfoIntegerArray(std::string(token.value().first));
      if (integer_array.size() != array_size) {

        ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Integer Array size : {} differs from token size: {}, field ID:{}, Number:{}, Type:{}",
                                integer_array.size(), array_size, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

      }

      // Add to the memory array.
      size_t memory_index = dynamic_accounting.integerCount();
      for (InfoIntegerType integer_value : integer_array) {

        integer_memory_[memory_index] = integer_value;
        ++memory_index;

      }

      // Update the dynamic accounting (this time for integers)
      dynamic_accounting.integerCountAdd(array_size);

    } else {

      integer_memory_[subscribed_field.dataOffset().offset] = VCFInfoParser::convertToInteger(std::string(token.value().first));

    }

  } else {

    integer_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_INTEGER;

  }

  return true;

}

bool kgl::InfoDataBlock::internUnityFloat(const InfoSubscribedField& subscribed_field,
                                          const std::optional<InfoParserToken>& token,
                                          InfoDataUsageCount& dynamic_accounting,
                                          InfoDataUsageCount& static_accounting) {

  if (subscribed_field.dataOffset().offset != static_accounting.floatCount()) {

    ExecEnv::log().critical(
    "InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}, field ID:{}, Number:{}, Type:{}",
    subscribed_field.dataOffset().offset, static_accounting.floatCount(), subscribed_field.infoRecord().ID,
    subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

  }
  // Increment the accounting.
  static_accounting.floatCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second > 1) {

      // Put a missing value in the static storage allocated.
      float_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_FLOAT;

      if (dynamic_accounting.unityArrayCount() > type_count_.unityArrayCount()) {

        ExecEnv::log().error("InfoDataBlock::indexAndVerify, Dynamic Unity Array Size: {} exceeds allocated Unity Array Size: {}",
                             dynamic_accounting.unityArrayCount(), type_count_.unityArrayCount());
        return false;
      }
      // Insert the array details into the unity array object.
      // Note the field index is not used.
      size_t array_size = token.value().second;
      unity_array_memory_[dynamic_accounting.unityArrayCount()].infoOffset(dynamic_accounting.floatCount());
      unity_array_memory_[dynamic_accounting.unityArrayCount()].infoSize(array_size);
      unity_array_memory_[dynamic_accounting.unityArrayCount()].infoVariableIndex(subscribed_field.fieldIndex());

      // Update the unity array count.
      dynamic_accounting.unityArrayCountAdd(1);

      // Copy in the float data.
      if (dynamic_accounting.floatCount() + array_size > type_count_.floatCount()) {

        ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Float Array size+offset :{} has exceeded Float vector size: {}, field ID:{}, Number:{}, Type:{}",
                                dynamic_accounting.floatCount() + array_size, type_count_.floatCount(), subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number,
                                subscribed_field.infoRecord().type);

      }
      // Get a vector of floats
      InfoParserFloatArray float_array =  VCFInfoParser::getInfoFloatArray(std::string(token.value().first));
      if (float_array.size() != array_size) {

        ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Float Array size : {} differs from token size: {}, field ID:{}, Number:{}, Type:{}",
                                float_array.size(), array_size, subscribed_field.infoRecord().ID, subscribed_field.infoRecord().number, subscribed_field.infoRecord().type);

      }

      // Add to the memory array.
      size_t memory_index = dynamic_accounting.floatCount();
      for (InfoFloatType float_value : float_array) {

        float_memory_[memory_index] = float_value;
        ++memory_index;

      }

      // Update the dynamic accounting (this time for integers)
      dynamic_accounting.floatCountAdd(array_size);

    } else {

      float_memory_[subscribed_field.dataOffset().offset] = VCFInfoParser::convertToFloat(std::string(token.value().first));

    }

  } else {

    float_memory_[subscribed_field.dataOffset().offset] = MISSING_VALUE_FLOAT;

  }

  return true;

}

