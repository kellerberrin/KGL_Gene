//
// Created by kellerberrin on 25/5/20.
//


#include "kgl_variant_factory_vcf_evidence_data_mem.h"



namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//  Returns the field size, A size of zero is a missing field.
size_t kgl::InfoDataBlock::getDataSize( size_t field_address,     // The field index
                                        size_t field_id,           // the field identifier
                                        InfoEvidenceIntern internal_type) const {


  switch (internal_type) {

    // Data is pre-allocated for the fixed fields.
    case InfoEvidenceIntern::intern_char:
      return 1;

    case InfoEvidenceIntern::intern_integer:
      return getInteger(field_address) != MISSING_VALUE_INTEGER_ ? 1 : 0;

    case InfoEvidenceIntern::intern_float:
      return getFloat(field_address) != MISSING_VALUE_FLOAT_ ? 1 : 0;

    case InfoEvidenceIntern::intern_unity_integer_array:
      return getIntegerUnityArray(field_address, field_id).infoSize();

    case InfoEvidenceIntern::intern_unity_float_array:
      return getFloatUnityArray(field_address, field_id).infoSize();

    case InfoEvidenceIntern::intern_string:
      return not getStringView(field_address).empty() ? 1 : 0;

    case InfoEvidenceIntern::intern_integer_array:
    case InfoEvidenceIntern::intern_float_array:
    case InfoEvidenceIntern::intern_unity_string_array:
    case InfoEvidenceIntern::intern_string_array:
      return getArray(field_address, field_id).infoSize();

    case InfoEvidenceIntern::NotImplemented:  // unknown internal type.
    default:
      ExecEnv::log().error( "InfoDataBlock::getDataSize, Internal data type unknown, cannot get size");
      return 0;

  }

}


bool kgl::InfoDataBlock::getBoolean(size_t field_address) const {

  if (field_address > type_count_.charCount()) {

    ExecEnv::log().critical("InfoDataBlock::getBoolean, Field Offset : {} exceeds char vector size: {}", field_address, type_count_.charCount());

  }

  return char_memory_[field_address];

}

kgl::InfoIntegerType kgl::InfoDataBlock::getInteger(size_t field_address) const {

  if (field_address > type_count_.integerCount()) {

    ExecEnv::log().critical("InfoDataBlock::getInteger, Field Offset : {} exceeds integer vector size: {}", field_address, type_count_.integerCount());

  }

  return integer_memory_[field_address];

}


kgl::InfoFloatType kgl::InfoDataBlock::getFloat(size_t field_address) const {

  if (field_address > type_count_.floatCount()) {

    ExecEnv::log().critical("InfoDataBlock::getFloat, Field Offset : {} exceeds float vector size: {}", field_address, type_count_.floatCount());

  }

  return float_memory_[field_address];

}


std::string_view kgl::InfoDataBlock::getStringView(size_t field_address) const {

  if (field_address > type_count_.stringCount()) {

    ExecEnv::log().critical("InfoDataBlock::getStringView, Field Offset : {} exceeds string vector size: {}", field_address, type_count_.stringCount());

  }

  return string_memory_[field_address];

}


kgl::InfoArrayIndex kgl::InfoDataBlock::getArray(size_t field_address, size_t field_id) const {

  if (field_address > type_count_.arrayCount()) {

    ExecEnv::log().critical("InfoDataBlock::getArray, Field Offset : {} exceeds array block size: {}", field_address, type_count_.arrayCount());

  }

  InfoArrayIndex array_block = array_memory_[field_address];

  if (array_block.infoVariableIndex() != field_id) {

    ExecEnv::log().critical("InfoDataBlock::getArray, Field Identifier : {} does not match array block field identifer: {}", field_id, array_block.infoVariableIndex());

  }

  return array_block;

}

// Unity arrays are used for Type='A' (alternate allele) variables.
// Most of these variables have a single value, but some (a few) are arrays with values for each alternate allele.
// Since parsimonious memory utilization was the primary consideration for the design of the Info data memory block.
// We initially assume that Type='A' variables will only use 1 element. If we examine this element and find a missing value
// Then we go to the unity_array_memory_ vector and look for an array block that matches the field_id. If this is not
// found them we assume that the variable is missing.
kgl::InfoArrayIndex kgl::InfoDataBlock::getIntegerUnityArray(size_t field_address, size_t field_id) const {

  InfoIntegerType unity_integer = getInteger(field_address);

  // Search for a match array block
  if (unity_integer == MISSING_VALUE_INTEGER_) {

    // search the unity array block
    for (size_t index = 0; index < type_count_.unityArrayCount(); ++index) {

      if (field_id == unity_array_memory_[index].infoVariableIndex()) {

        return unity_array_memory_[index];

      }

    }

    // Debug.
    // ExecEnv::log().warn("InfoDataBlock::getIntegerUnityArray, (Debug Remove) Missing value Field id: {}, Field Address: {}", field_id, field_address);
    InfoArrayIndex empty_sized_array;
    empty_sized_array.infoVariableIndex(field_id);
    empty_sized_array.infoOffset(0);
    empty_sized_array.infoSize(0);

    // return empty array block.
    return empty_sized_array;

  } else {

    // Construct an array block size 1 and return.
    InfoArrayIndex unit_sized_array;
    unit_sized_array.infoVariableIndex(field_id);
    unit_sized_array.infoOffset(field_address);
    unit_sized_array.infoSize(1);

    return unit_sized_array;

  }

}

kgl::InfoArrayIndex kgl::InfoDataBlock::getFloatUnityArray(size_t field_address, size_t field_id) const {

  InfoFloatType unity_float = getFloat(field_address);

  // Search for a match array block
  if (unity_float == MISSING_VALUE_FLOAT_) {

    // search the unity array block
    for (size_t index = 0; index < type_count_.unityArrayCount(); ++index) {

      if (field_id == unity_array_memory_[index].infoVariableIndex()) {

        return unity_array_memory_[index];

      }

    }

    // Debug.
    // ExecEnv::log().warn("InfoDataBlock::getFloatUnityArray, (Debug Remove) Missing value Field id: {}, Field Address: {}", field_id, field_address);
    InfoArrayIndex empty_sized_array;
    empty_sized_array.infoVariableIndex(field_id);
    empty_sized_array.infoOffset(0);
    empty_sized_array.infoSize(0);

    // return empty array block.
    return empty_sized_array;

  } else {

    // Construct an array block size 1 and return.
    InfoArrayIndex unit_sized_array;
    unit_sized_array.infoVariableIndex(field_id);
    unit_sized_array.infoOffset(field_address);
    unit_sized_array.infoSize(1);

    return unit_sized_array;

  }

}

std::string kgl::InfoDataBlock::getString(size_t field_address) const {

  std::string_view view = getStringView(field_address);

  auto char_base = reinterpret_cast<size_t>(&char_memory_[0]);

  auto char_offset = reinterpret_cast<size_t>(view.data());

  if (not view.empty() and (char_offset-char_base) + view.size() > type_count_.charCount()) {

    ExecEnv::log().critical("InfoDataBlock::getString, string Offset: {} + String Size: {}, exceeds the char vector size: {}",
                            (char_offset-char_base), view.size(), type_count_.charCount());

  }

  return std::string(view);

}


std::vector<int64_t> kgl::InfoDataBlock::getIntegerArray(size_t field_address, size_t field_id) const {

  InfoArrayIndex array_index = getArray(field_address, field_id);

  if (array_index.infoOffset() + array_index.infoSize() > type_count_.integerCount()) {

    ExecEnv::log().critical( "InfoDataBlock::getIntegerArray, Array Offset: {} + Array Size: {}  exceeds exceeds integer vector size: {}",
                             array_index.infoOffset(), array_index.infoSize(), type_count_.integerCount());

  }
  std::vector<int64_t> integer_vector;
  integer_vector.reserve(array_index.infoSize());
  for (size_t index = 0; index < array_index.infoSize(); ++index) {

    integer_vector.push_back(static_cast<int64_t>(integer_memory_[(array_index.infoOffset() + index)]));

  }

  return integer_vector;

}


std::vector<double> kgl::InfoDataBlock::getFloatArray(size_t field_address, size_t field_id) const {

  InfoArrayIndex array_index = getArray(field_address, field_id);

  if (array_index.infoOffset() + array_index.infoSize() > type_count_.floatCount()) {

    ExecEnv::log().critical( "InfoDataBlock::getFloatArray, Array Offset: {} + Array Size: {}  exceeds exceeds float vector size: {}",
                             array_index.infoOffset(), array_index.infoSize(), type_count_.floatCount());

  }
  std::vector<double> float_vector;
  float_vector.reserve(array_index.infoSize());
  for (size_t index = 0; index < array_index.infoSize(); ++index) {

    float_vector.push_back(static_cast<double>(float_memory_[(array_index.infoOffset() + index)]));

  }

  return float_vector;

}


std::vector<std::string> kgl::InfoDataBlock::getStringArray(size_t field_address, size_t field_id) const {


  InfoArrayIndex array_index = getArray(field_address, field_id);

  if (array_index.infoOffset() + array_index.infoSize() > type_count_.stringCount()) {

    ExecEnv::log().critical( "InfoDataBlock::getStringArray, Array Offset: {} + Array Size: {}  exceeds exceeds string vector size: {}",
                             array_index.infoOffset(), array_index.infoSize(), type_count_.stringCount());

  }
  std::vector<std::string_view> string_view_vector;
  string_view_vector.reserve(array_index.infoSize());
  for (size_t index = 0; index < array_index.infoSize(); ++index) {

    string_view_vector.push_back(string_memory_[(array_index.infoOffset() + index)]);

  }

  std::vector<std::string> string_vector;
  string_vector.reserve(array_index.infoSize());
  auto char_base = reinterpret_cast<size_t>(&char_memory_[0]);
  for (auto const& view : string_view_vector) {

    auto char_offset = reinterpret_cast<size_t>(view.data());

    if (not view.empty() and (char_offset-char_base) + view.size() > type_count_.charCount()) {

      ExecEnv::log().critical("InfoDataBlock::getStringArray, string Offset: {} + String Size: {}, exceeds the char vector size: {}",
                              (char_offset-char_base), view.size(), type_count_.charCount());

    }

    string_vector.emplace_back(std::string(view));

  }

  return string_vector;

}


std::vector<int64_t> kgl::InfoDataBlock::getUnityIntegerArray(size_t field_address, size_t field_id) const {

  InfoArrayIndex array_index = getIntegerUnityArray(field_address, field_id);

  if (array_index.infoOffset() + array_index.infoSize() > type_count_.integerCount()) {

    ExecEnv::log().critical( "InfoDataBlock::getUnityIntegerArray, Array Offset: {} + Array Size: {}  exceeds exceeds integer vector size: {}",
                             array_index.infoOffset(), array_index.infoSize(), type_count_.integerCount());

  }
  std::vector<int64_t> integer_vector;
  integer_vector.reserve(array_index.infoSize());
  for (size_t index = 0; index < array_index.infoSize(); ++index) {

    integer_vector.push_back(static_cast<int64_t>(integer_memory_[(array_index.infoOffset() + index)]));

  }

  return integer_vector;

}


std::vector<double> kgl::InfoDataBlock::getUnityFloatArray(size_t field_address, size_t field_id) const {


  InfoArrayIndex array_index = getFloatUnityArray(field_address, field_id);

  if (array_index.infoOffset() + array_index.infoSize() > type_count_.floatCount()) {

    ExecEnv::log().critical( "InfoDataBlock::getUnityFloatArray, Array Offset: {} + Array Size: {}  exceeds exceeds float vector size: {}",
                             array_index.infoOffset(), array_index.infoSize(), type_count_.floatCount());

  }
  std::vector<double> float_vector;
  float_vector.reserve(array_index.infoSize());
  for (size_t index = 0; index < array_index.infoSize(); ++index) {

    float_vector.push_back(static_cast<double>(float_memory_[(array_index.infoOffset() + index)]));

  }

  return float_vector;

}




