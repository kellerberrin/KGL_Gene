//
// Created by kellerberrin on 24/5/20.
//

#include "kgl_variant_factory_vcf_evidence_data_mem.h"


namespace kgl = kellerberrin::genome;


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


// Actually stores the Info data in the raw memory and does index and usage checking.
bool kgl::InfoDataBlock::indexAndVerify( size_t field_address,     // The field index
                                         size_t field_id,           // the field identifier
                                         InfoEvidenceIntern internal_type,  // The internal data type.
                                         const std::optional<InfoParserToken>& token, // The parsed token to be placed in memory
                                         InfoDataUsageCount& dynamic_accounting,  // accounting objects to check index, data and memory integrity
                                         InfoDataUsageCount& static_accounting) {


  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      return internChar(field_address, token, static_accounting);

    case InfoEvidenceIntern::intern_integer:
      return internInteger(field_address, token, static_accounting);

    case InfoEvidenceIntern::intern_float:
      return internFloat(field_address, token, static_accounting);

    case InfoEvidenceIntern::intern_unity_integer_array:
      return internUnityInteger(field_address, field_id, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_unity_float_array:
      return internUnityFloat(field_address, field_id, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_string:
      return internString(field_address, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_integer_array:
      return internIntegerArray(field_address, field_id, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_float_array:
      return internFloatArray(field_address, field_id, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array:
      return internStringArray(field_address, field_id, token, dynamic_accounting, static_accounting);

    case InfoEvidenceIntern::NotImplemented:  // unknown internal type.
    default:
      ExecEnv::log().error( "InfoDataBlock::indexAndVerify, Internal data type unknown");
      return false;

  }

}


// These are just function blocks that break up the unsightly mass of code in function InfoDataBlock::indexAndVerify.
bool kgl::InfoDataBlock::internChar(size_t field_address,
                                    const std::optional<InfoParserToken>& token, // The parsed token to be placed in memory
                                    InfoDataUsageCount& static_accounting) {

  // check field index
  if (field_address != static_accounting.charCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.charCount());
    return false;

  }
  // Increment accounting.
  static_accounting.charCountAdd(1);
  // Add the data.
  char_memory_[field_address] = token ? true : false;

  return true;

}

bool kgl::InfoDataBlock::internInteger(size_t field_address,
                                       const std::optional<InfoParserToken>& token,
                                       InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.integerCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.integerCount());
    return false;

  }

  // Increment the accounting.
  static_accounting.integerCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}", token.value().second);
      return false;

    }

    integer_memory_[field_address] = VCFInfoParser::convertToInteger(std::string(token.value().first));

  } else {

    integer_memory_[field_address] = MISSING_VALUE_INTEGER_;

  }

  return true;

}

bool kgl::InfoDataBlock::internFloat(size_t field_address,
                                     const std::optional<InfoParserToken>& token,
                                     InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.floatCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.floatCount());
    return false;

  }
  // Increment the accounting.
  static_accounting.floatCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, Scalar Field has vector size: {}", token.value().second);
      return false;

    }

    float_memory_[field_address] = VCFInfoParser::convertToFloat(std::string(token.value().first));

  } else {

    float_memory_[field_address] = MISSING_VALUE_FLOAT_;

  }

  return true;

}

bool kgl::InfoDataBlock::internString(size_t field_address,
                                      const std::optional<InfoParserToken>& token,
                                      InfoDataUsageCount& dynamic_accounting,
                                      InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.stringCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.stringCount());
    return false;

  }
  // Update the static accounting
  static_accounting.stringCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second != 1) {

      // A scalar string can have embedded array delimiter "," characters which confuses the parser.
      ExecEnv::log().warn("InfoDataBlock::indexAndVerify, String field has parsed size: {} (expected 1) contains ',' special characters", token.value().second);

    }

    char* string_start = &char_memory_[dynamic_accounting.charCount()];
    size_t string_size = token.value().first.size();
    string_memory_[field_address] = std::string_view(string_start, string_size);
    if (dynamic_accounting.charCount() + string_size > type_count_.charCount()) {

      ExecEnv::log().error( "InfoDataBlock::indexAndVerify, String size+offset :{} has exceeded char vector size: {}",
                            dynamic_accounting.charCount() + string_size, type_count_.charCount());
      return false;

    }
    std::memcpy(string_start, token.value().first.data(), string_size);
    // Update the dynamic accounting.
    dynamic_accounting.charCountAdd(string_size);

  } else {
// Empty string.
    string_memory_[field_address] = std::string_view();

  }

  return true;

}


bool kgl::InfoDataBlock::internIntegerArray(size_t field_address,
                                            size_t field_id,
                                            const std::optional<InfoParserToken>& token,
                                            InfoDataUsageCount& dynamic_accounting,
                                            InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.arrayCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.arrayCount());
    return false;

  }
  // Update the static accounting
  static_accounting.arrayCountAdd(1);
  // add the data.
  if (token) {

    size_t array_size = token.value().second;
    array_memory_[field_address] = InfoArrayIndex(field_id, dynamic_accounting.integerCount(), array_size);

    if (dynamic_accounting.integerCount() + array_size > type_count_.integerCount()) {

      ExecEnv::log().error( "InfoDataBlock::indexAndVerify, Integer Array offset+size :{} has exceeded Integer vector size: {}",
                            dynamic_accounting.integerCount() + array_size, type_count_.integerCount());
      return false;
    }
    // Get a vector of integers
    InfoParserIntegerArray integer_array =  VCFInfoParser::getInfoIntegerArray(std::string(token.value().first));
    if (integer_array.size() != array_size) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, Integer Array size : {} differs from token size: {}", integer_array.size(), array_size);
      return false;

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
    array_memory_[field_address] = InfoArrayIndex(field_id, 0, 0);

  }

  return true;

}


bool kgl::InfoDataBlock::internFloatArray(size_t field_address,
                                          size_t field_id,
                                          const std::optional<InfoParserToken>& token,
                                          InfoDataUsageCount& dynamic_accounting,
                                          InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.arrayCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}",
                         field_address, static_accounting.arrayCount());
    return false;
  }
  // Update the static accounting
  static_accounting.arrayCountAdd(1);
  // add the data.
  if (token) {

    size_t array_size = token.value().second;
    array_memory_[field_address] = InfoArrayIndex(field_id, dynamic_accounting.floatCount(), array_size);

    if (dynamic_accounting.floatCount() + array_size > type_count_.floatCount()) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, Float Array size+offset :{} has exceeded Float vector size: {}",
                           dynamic_accounting.floatCount() + array_size, type_count_.floatCount());
      return false;

    }
    // Get a vector of integers
    InfoParserFloatArray float_array =  VCFInfoParser::getInfoFloatArray(std::string(token.value().first));
    if (float_array.size() != array_size) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, Float Array size : {} differs from token size: {}", float_array.size(), array_size);
      return false;

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
    array_memory_[field_address] = InfoArrayIndex(field_id, 0, 0);

  }

  return true;

}


bool kgl::InfoDataBlock::internStringArray(size_t field_address,
                                           size_t field_id,
                                           const std::optional<InfoParserToken>& token,
                                           InfoDataUsageCount& dynamic_accounting,
                                           InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.arrayCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.arrayCount());
    return false;

  }
  // Update the static accounting
  static_accounting.arrayCountAdd(1);
  // add the data.
  if (token) {

    size_t array_size = token.value().second;
    array_memory_[field_address] = InfoArrayIndex(field_id, dynamic_accounting.stringCount(), array_size);

    if (dynamic_accounting.stringCount() + array_size > type_count_.stringCount()) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, String Array size+offset :{} has exceeded String vector size: {}",
                           dynamic_accounting.stringCount() + array_size, type_count_.stringCount());
      return false;

    }
    // Get a vector of string_views.
    std::vector<std::string_view> string_view_array = VCFInfoParser::getInfoStringArray(token.value().first);
    if (string_view_array.size() != array_size) {

      ExecEnv::log().error("InfoDataBlock::indexAndVerify, String Array size : {} differs from token size: {}", string_view_array.size(), array_size);
      return false;

    }

    // Add to the memory array.
    size_t memory_index = dynamic_accounting.stringCount();
    for (auto str_view : string_view_array) {

      char* string_address = &char_memory_[dynamic_accounting.charCount()];
      string_memory_[memory_index] = std::string_view(string_address, str_view.size());
      // update char count.
      dynamic_accounting.charCountAdd(str_view.size());

      if (dynamic_accounting.charCount() > type_count_.charCount()) {

        ExecEnv::log().error( "InfoDataBlock::indexAndVerify, String size+offset :{} has exceeded char vector size: {}",
                              dynamic_accounting.charCount(), type_count_.charCount());
        return false;

      }

      std::memcpy(string_address, str_view.data(), str_view.size());
      ++memory_index;

    }

    // Update the dynamic accounting for strings.
    dynamic_accounting.stringCountAdd(array_size);

  } else {
// Empty vector
    array_memory_[field_address] = InfoArrayIndex(field_id, 0, 0);

  }

  return true;

}

bool kgl::InfoDataBlock::internUnityInteger(size_t field_address,
                                            size_t field_id,
                                            const std::optional<InfoParserToken>& token,
                                            InfoDataUsageCount& dynamic_accounting,
                                            InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.integerCount()) {

    ExecEnv::log().error("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}", field_address, static_accounting.integerCount());
    return false;

  }

  // Increment the static accounting.
  static_accounting.integerCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second > 1) {
// We have a vector, which can (infrequently) occur with this internal datatype.

      // Put a missing value in the static storage allocated.
      integer_memory_[field_address] = MISSING_VALUE_INTEGER_;

      if (dynamic_accounting.unityArrayCount() > type_count_.unityArrayCount()) {

        ExecEnv::log().error("InfoDataBlock::indexAndVerify, Dynamic Unity Array Size: {} exceeds allocated Unity Array Size: {}",
                             dynamic_accounting.unityArrayCount(), type_count_.unityArrayCount());
        return false;
      }
      // Insert the array details into the unity array object.
      // Note the field index is not used.
      size_t array_size = token.value().second;
      unity_array_memory_[dynamic_accounting.unityArrayCount()] = InfoArrayIndex(field_id, dynamic_accounting.integerCount(), array_size);

      // Update the unity array count.
      dynamic_accounting.unityArrayCountAdd(1);

      // Copy in the integer data.
      if (dynamic_accounting.integerCount() + array_size > type_count_.integerCount()) {

        ExecEnv::log().error("InfoDataBlock::indexAndVerify, Integer Array size+offset :{} has exceeded Integer vector size: {}",
                             dynamic_accounting.integerCount() + array_size, type_count_.integerCount());
        return false;

      }
      // Get a vector of integers
      InfoParserIntegerArray integer_array =  VCFInfoParser::getInfoIntegerArray(std::string(token.value().first));
      if (integer_array.size() != array_size) {

        ExecEnv::log().error("InfoDataBlock::indexAndVerify, Integer Array size : {} differs from token size: {}", integer_array.size(), array_size);
        return false;

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

      // A scalar integer
      integer_memory_[field_address] = VCFInfoParser::convertToInteger(std::string(token.value().first));

    }

  } else {
    // Actual missing value.
    integer_memory_[field_address] = MISSING_VALUE_INTEGER_;

  }

  return true;

}

bool kgl::InfoDataBlock::internUnityFloat(size_t field_address,
                                          size_t field_id,
                                          const std::optional<InfoParserToken>& token,
                                          InfoDataUsageCount& dynamic_accounting,
                                          InfoDataUsageCount& static_accounting) {

  if (field_address != static_accounting.floatCount()) {

    ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Field offset: {} differs from calculated offset: {}",
                            field_address, static_accounting.floatCount());

  }
  // Increment the accounting.
  static_accounting.floatCountAdd(1);
  // add the data.
  if (token) {

    if (token.value().second > 1) {

      // Put a missing value in the static storage allocated.
      float_memory_[field_address] = MISSING_VALUE_FLOAT_;

      if (dynamic_accounting.unityArrayCount() > type_count_.unityArrayCount()) {

        ExecEnv::log().error("InfoDataBlock::indexAndVerify, Dynamic Unity Array Size: {} exceeds allocated Unity Array Size: {}",
                             dynamic_accounting.unityArrayCount(), type_count_.unityArrayCount());
        return false;
      }
      // Insert the array details into the unity array object.
      // Note the field index is not used.
      size_t array_size = token.value().second;
      unity_array_memory_[dynamic_accounting.unityArrayCount()] = InfoArrayIndex(field_id, dynamic_accounting.floatCount(), array_size);

      // Update the unity array count.
      dynamic_accounting.unityArrayCountAdd(1);

      // Copy in the float data.
      if (dynamic_accounting.floatCount() + array_size > type_count_.floatCount()) {

        ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Float Array size+offset :{} has exceeded Float vector size: {}",
                                dynamic_accounting.floatCount() + array_size, type_count_.floatCount());

      }
      // Get a vector of floats
      InfoParserFloatArray float_array =  VCFInfoParser::getInfoFloatArray(std::string(token.value().first));
      if (float_array.size() != array_size) {

        ExecEnv::log().critical("InfoDataBlock::indexAndVerify, Float Array size : {} differs from token size: {}",
                                float_array.size(), array_size);

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

      float_memory_[field_address] = VCFInfoParser::convertToFloat(std::string(token.value().first));

    }

  } else {

    float_memory_[field_address] = MISSING_VALUE_FLOAT_;

  }

  return true;

}

