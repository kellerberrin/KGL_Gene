//
// Created by kellerberrin on 30/5/20.
//

#include "kgl_variant_factory_vcf_evidence_data_blk.h"
#include "kgl_variant_factory_vcf_evidence.h"

#include <algorithm>

namespace kgl = kellerberrin::genome;


kgl::DataMemoryBlock::DataMemoryBlock( std::shared_ptr<const InfoEvidenceHeader> info_evidence_header,
                                       const InfoMemoryResource& memory_resource,
                                       const VCFInfoParser& info_parser)
                                       : info_evidence_header_(std::move(info_evidence_header)) {

  // booleans and strings are both implemented as type char
  type_count_.charCount(memory_resource.boolSize() + memory_resource.charSize());
  type_count_.floatCount(memory_resource.floatSize());
  type_count_.stringCount(memory_resource.viewSize());
  type_count_.integerCount(memory_resource.integerSize());
  type_count_.arrayCount(memory_resource.arraySize());
  // Perform actual memory allocation.
  char_memory_ = std::make_unique<char[]>(type_count_.charCount());
  integer_memory_ = std::make_unique<InfoIntegerType[]>(type_count_.integerCount());
  float_memory_ = std::make_unique<InfoFloatType[]>(type_count_.floatCount()) ;
  array_memory_ = std::make_unique<InfoArrayIndex[]>(type_count_.arrayCount());
  string_memory_ = std::make_unique<std::string_view[]>(type_count_.stringCount());
  // Copy the array blocks.
  // The dynamicAllocation() map ensures that copied array memory blocks will be sorted in ascending identifier order.
  // We can use a binary search for fast access.
  size_t index = 0;
  for (auto const& array_block : memory_resource.arrayResource().dynamicAllocation()) {

    array_memory_[index] = array_block.second;
    ++index;

  }
  // Resource instance keeps track of string char usage.
  FixedResourceInstance string_usage("Char Counter");
  // allocate the first part of the char array to booleans.
  string_usage.allocateResource(memory_resource.boolSize());
  // Store all the parser data in the memory block.
  for (auto const& [data_identifier, data_item] : info_evidence_header_->getConstMap()) {

    // Lookup the matching data token.
    std::optional<const InfoParserToken> token = info_parser.getToken(data_identifier);

    switch(data_item.getDataHandle()->resourceType()) {

      case DataResourceType::Boolean:
        storeBoolean(memory_resource, data_item.getDataHandle().value(), token);
        break;

      case DataResourceType::Integer:
        storeInteger(data_item.getDataHandle().value(), token);
        break;

      case DataResourceType::Float:
        storeFloat(data_item.getDataHandle().value(), token);
        break;

      case DataResourceType::String:
        storeString(data_item.getDataHandle().value(), token, string_usage);
        break;

    }

  }

}

// Binary lookup on the array index.
const kgl::InfoArrayIndex& kgl::DataMemoryBlock::findIdentifier(size_t identifier) const {

  const auto comp = [](const InfoArrayIndex& a, const InfoArrayIndex& b) -> bool { return a.infoVariableIndex() < b.infoVariableIndex(); };

  auto result_pair = std::equal_range(&array_memory_[0],&array_memory_[type_count_.arrayCount()], InfoArrayIndex(identifier, 0, 0), comp);

  if (result_pair.first->infoVariableIndex() == identifier) {

    return *result_pair.first;

  } else {

    ExecEnv::log().error("DataMemoryBlock::findIdentifier, Identifier lookup failed");
    return *result_pair.first;

  }

}



void kgl::DataMemoryBlock::storeBoolean(const InfoMemoryResource& memory_resource, const InfoResourceHandle& handle, std::optional<const InfoParserToken> token) {

  // Check that the bool is fixed with size = 1.
  if (handle.dynamicType() != DataDynamicType::FixedData or handle.initialDataSize() != 1) {

    ExecEnv::log().error("DataMemoryBlock::storeBoolean, Boolean type should be Fixed size = 1");
    return;

  }

  // Check that the offset does not exceed the bounds of the boolean block.

  if (handle.initialDataOffset() >= memory_resource.boolSize()) {

    ExecEnv::log().error( "DataMemoryBlock::storeBoolean, Data Offset: {} exceeds maximum calculated boolean offset: {}",
                          handle.initialDataOffset(), memory_resource.boolSize());
    return;

  }

  char_memory_[handle.initialDataOffset()] = token ? true : false;

}

void kgl::DataMemoryBlock::storeString(const InfoResourceHandle& handle, std::optional<const InfoParserToken> token, FixedResourceInstance& string_usage) {

  if (handle.dynamicType() == DataDynamicType::FixedData) {

    if (handle.initialDataOffset() >= type_count_.stringCount()) {

      ExecEnv::log().error( "DataMemoryBlock::storeString, Data Offset: {} exceeds maximum calculated string view offset: {}",
                            handle.initialDataOffset(), type_count_.stringCount());
      return;

    }

    // No token means a string size of zero (empty string)
    size_t string_size = token ? token.value().first.size() : 0;
    size_t string_offset = string_usage.allocateResource(string_size); // increments to the next string.

    if (string_usage.resourceValue() > type_count_.charCount()) {

      ExecEnv::log().error( "DataMemoryBlock::storeString, Data Offset: {} exceeds maximum calculated string view offset: {}",
                            string_usage.resourceValue(), type_count_.charCount());
      return;

    }

    // Get the destination address
    char* string_dest_ptr = &char_memory_[string_offset];
    // store the string_view
    string_memory_[handle.initialDataOffset()] = std::string_view(string_dest_ptr, string_size);
    // If valid token, then copy the string.
    if (token) {

      std::memcpy(string_dest_ptr, token.value().first.data(), string_size);

    }

  } else if (handle.dynamicType() == DataDynamicType::DynamicData) {

    const InfoArrayIndex& array_index = findIdentifier(handle.handleId());

    // Ensure we have space in the string view vector for all the strings.
    if (array_index.infoOffset() + array_index.infoSize() > type_count_.stringCount()) {

      ExecEnv::log().error(
      "DataMemoryBlock::storeString, Identifier: {} Array Offset: {} + Size: {} ({})exceeds maximum calculated string view offset: {}",
      handle.handleId(), array_index.infoOffset(), array_index.infoSize(), (array_index.infoOffset() + array_index.infoSize()), type_count_.stringCount());
      return;

    }

    // Check if we have a valid token.
    if (token) {
      // Check the token size and the array size.
      if (array_index.infoSize() != token.value().second) {

        ExecEnv::log().error("DataMemoryBlock::storeString, Mismatch between array size {} and token size: {}",
                             array_index.infoSize(), array_index.infoSize());
        return;

      }

      // Parse out the strings.
      std::vector<std::string_view> string_vector = Utility::view_tokenizer(token.value().first,
                                                                            VCFInfoParser::INFO_VECTOR_DELIMITER_);

      // Ensure the array block is the same size as the parsed string count.
      if (array_index.infoSize() != string_vector.size()) {

        ExecEnv::log().error("DataMemoryBlock::storeString, Array Size: {} not equal to string array size: {}",
                             array_index.infoSize(), string_vector.size());
        return;

      }

      // Copy each string.
      size_t index = array_index.infoOffset();
      for (auto const &str_view : string_vector) {

        size_t string_offset = string_usage.allocateResource(str_view.size()); // increments to the next string.

        if (string_usage.resourceValue() > type_count_.charCount()) {

          ExecEnv::log().error(
          "DataMemoryBlock::storeString, Data Offset: {} exceeds maximum calculated string view offset: {}",
          string_usage.resourceValue(), type_count_.charCount());
          return;

        }

        // Get the destination address
        char *string_dest_ptr = &char_memory_[string_offset];
        // store the string_view
        string_memory_[index] = std::string_view(string_dest_ptr, str_view.size());
        // copy the string.
        std::memcpy(string_dest_ptr, str_view.data(), str_view.size());
        // next view
        ++index;

      }

    } else { //token not valid so assign zero-sized strings to all array elements.

      for (size_t index = 0; index < array_index.infoSize(); ++index) {

        string_memory_[(index + array_index.infoOffset())] = std::string_view(nullptr, 0);

      }

    }

  } else {

    ExecEnv::log().error( "DataMemoryBlock::storeString, Strings only supported as fixed or dynamic data types");
    return;

  }

}

void kgl::DataMemoryBlock::storeInteger(const InfoResourceHandle& handle, std::optional<const InfoParserToken> token) {

  if (handle.dynamicType() == DataDynamicType::FixedData and handle.initialDataSize() == 1) {

    if (handle.initialDataOffset() >= type_count_.integerCount()) {

      ExecEnv::log().error( "DataMemoryBlock::storeInteger, Data Offset: {} exceeds maximum calculated integer offset: {}",
                            handle.initialDataOffset(), type_count_.integerCount());
      return;

    }

    if (token) {

      if (token.value().second > 1) {

        ExecEnv::log().error( "DataMemoryBlock::storeInteger, Expected data size 1, token size: {}, token value: {}",
                              token.value().second, std::string(token.value().first));
        return;

      }

      integer_memory_[handle.initialDataOffset()] = VCFInfoParser::convertToInteger(std::string(token.value().first));

    } else { // Token is not valid.

      integer_memory_[handle.initialDataOffset()] = MISSING_VALUE_INTEGER_;

    }

  } if (handle.dynamicType() == DataDynamicType::FixedDynamic) {

    if (handle.initialDataOffset() >= type_count_.integerCount()) {

      ExecEnv::log().error( "DataMemoryBlock::storeInteger, Data Offset: {} exceeds maximum calculated integer offset: {}",
                            handle.initialDataOffset(), type_count_.integerCount());
      return;

    }

    if (token) {

      if (token.value().second > 1) {  // Implemented as an array.

        const InfoArrayIndex& array_index = findIdentifier(handle.handleId());

        if (array_index.infoOffset() + array_index.infoSize() > type_count_.integerCount()) {

          ExecEnv::log().error("DataMemoryBlock::storeInteger, Array Offset: {} + Array Size: {} exceeds Integer Array Size: {}",
                               array_index.infoOffset(), array_index.infoSize(), type_count_.integerCount());
          return;

        }

        // Parse out the strings.
        std::vector<std::string_view> string_vector = Utility::view_tokenizer(token.value().first,
                                                                              VCFInfoParser::INFO_VECTOR_DELIMITER_);

        if (array_index.infoSize() != string_vector.size()) {

          ExecEnv::log().error("DataMemoryBlock::storeInteger, Mismatch between array size {} and token size: {}",
                               array_index.infoSize(), string_vector.size());
          return;

        }

        // Invalidate scalar
        integer_memory_[handle.initialDataOffset()] = MISSING_VALUE_INTEGER_;
        // Populate array
        size_t index = array_index.infoOffset();
        for (auto const& str_view : string_vector) {

          integer_memory_[index] = VCFInfoParser::convertToInteger(std::string(str_view));

          ++index;
        }

      } else { // Implemented as a scalar.

        integer_memory_[handle.initialDataOffset()] = VCFInfoParser::convertToInteger(std::string(token.value().first));

      }

    } else { // Token not valid.

      integer_memory_[handle.initialDataOffset()] = MISSING_VALUE_INTEGER_;

    }

  } if (handle.dynamicType() == DataDynamicType::DynamicData) {

    const InfoArrayIndex &array_index = findIdentifier(handle.handleId());

    if (array_index.infoOffset() + array_index.infoSize() > type_count_.integerCount()) {

      ExecEnv::log().error(
      "DataMemoryBlock::storeInteger, Array Offset: {} + Array Size: {} exceeds Integer Array Size: {}",
      array_index.infoOffset(), array_index.infoSize(), type_count_.integerCount());
      return;

    }

    if (token) {

      // Parse out the strings.
      std::vector<std::string_view> string_vector = Utility::view_tokenizer(token.value().first,
                                                                            VCFInfoParser::INFO_VECTOR_DELIMITER_);

      if (array_index.infoSize() != string_vector.size()) {

        ExecEnv::log().error("DataMemoryBlock::storeInteger, Mismatch between array size {} and token size: {}",
                             array_index.infoSize(), string_vector.size());
        return;

      }
      // Populate array.
      size_t index = array_index.infoOffset();
      for (auto const &str_view : string_vector) {

        integer_memory_[index] = VCFInfoParser::convertToInteger(std::string(str_view));

        ++index;
      }

    } else { // Token not valid.

      for (size_t index = 0; index < array_index.infoSize(); ++index) {

        integer_memory_[array_index.infoOffset() + index] = MISSING_VALUE_INTEGER_;

      } // for loop

    } // no token

  } // if dynamic

}

void kgl::DataMemoryBlock::storeFloat(const InfoResourceHandle& handle, std::optional<const InfoParserToken> token) {

  if (handle.dynamicType() == DataDynamicType::FixedData and handle.initialDataSize() == 1) {

    if (handle.initialDataOffset() >= type_count_.floatCount()) {

      ExecEnv::log().error( "DataMemoryBlock::storeFloat, Data Offset: {} exceeds maximum calculated integer offset: {}",
                            handle.initialDataOffset(), type_count_.floatCount());
      return;

    }

    if (token) {

      if (token.value().second > 1) {

        ExecEnv::log().error( "DataMemoryBlock::storeFloat, Expected data size 1, token size: {}, token value: {}",
                              token.value().second, std::string(token.value().first));
        return;

      }

      float_memory_[handle.initialDataOffset()] = VCFInfoParser::convertToFloat(std::string(token.value().first));

    } else { // Token is not valid.

      float_memory_[handle.initialDataOffset()] = MISSING_VALUE_FLOAT_;

    }

  } if (handle.dynamicType() == DataDynamicType::FixedDynamic) {

    if (handle.initialDataOffset() >= type_count_.floatCount()) {

      ExecEnv::log().error( "DataMemoryBlock::storeFloat, Data Offset: {} exceeds maximum calculated float offset: {}",
                            handle.initialDataOffset(), type_count_.floatCount());
      return;

    }

    if (token) {

      if (token.value().second > 1) {  // Implemented as an array.

        const InfoArrayIndex& array_index = findIdentifier(handle.handleId());

        if (array_index.infoOffset() + array_index.infoSize() > type_count_.floatCount()) {

          ExecEnv::log().error("DataMemoryBlock::storeFloat, Array Offset: {} + Array Size: {} exceeds Float Array Size: {}",
                               array_index.infoOffset(), array_index.infoSize(), type_count_.floatCount());
          return;

        }

        // Parse out the strings.
        std::vector<std::string_view> string_vector = Utility::view_tokenizer(token.value().first,
                                                                              VCFInfoParser::INFO_VECTOR_DELIMITER_);

        if (array_index.infoSize() != string_vector.size()) {

          ExecEnv::log().error("DataMemoryBlock::storeFloat, Mismatch between array size {} and token size: {}",
                               array_index.infoSize(), string_vector.size());
          return;

        }

        // Invalidate scalar.
        float_memory_[handle.initialDataOffset()] = MISSING_VALUE_FLOAT_;
        // Populate array
        size_t index = array_index.infoOffset();
        for (auto const& str_view : string_vector) {

          float_memory_[index] = VCFInfoParser::convertToFloat(std::string(str_view));

          ++index;
        }

      } else { // Implemented as a scalar.

        float_memory_[handle.initialDataOffset()] = VCFInfoParser::convertToFloat(std::string(token.value().first));

      }

    } else { // Token not valid.

      float_memory_[handle.initialDataOffset()] = MISSING_VALUE_FLOAT_;

    }

  } if (handle.dynamicType() == DataDynamicType::DynamicData) {

    const InfoArrayIndex &array_index = findIdentifier(handle.handleId());

    if (array_index.infoOffset() + array_index.infoSize() > type_count_.floatCount()) {

      ExecEnv::log().error(
      "DataMemoryBlock::storeFloat, Array Offset: {} + Array Size: {} exceeds Float Array Size: {}",
      array_index.infoOffset(), array_index.infoSize(), type_count_.floatCount());
      return;

    }

    if (token) {

      // Parse out the strings.
      std::vector<std::string_view> string_vector = Utility::view_tokenizer(token.value().first,
                                                                            VCFInfoParser::INFO_VECTOR_DELIMITER_);

      if (array_index.infoSize() != string_vector.size()) {

        ExecEnv::log().error("DataMemoryBlock::storeFloat, Mismatch between array size {} and token size: {}",
                             array_index.infoSize(), string_vector.size());
        return;

      }
      // Populate the array.
      size_t index = array_index.infoOffset();
      for (auto const &str_view : string_vector) {

        float_memory_[index] = VCFInfoParser::convertToFloat(std::string(str_view));

        ++index;
      }

    } else { // Token not valid.

      for (size_t index = 0; index < array_index.infoSize(); ++index) {

        float_memory_[array_index.infoOffset() + index] = MISSING_VALUE_FLOAT_;

      } // for loop

    } // no token

  } // if dynamic

}


