//
// Created by kellerberrin on 1/6/20.
//

#include "kgl_variant_factory_vcf_evidence_data_blk.h"


namespace kgl = kellerberrin::genome;


bool kgl::DataMemoryBlock::getBoolean(const InfoResourceHandle& handle) const {

  if (handle.initialDataOffset() >= mem_count_.charCount()) {

    ExecEnv::log().error("DataMemoryBlock::getBoolean, Handle offset: {} exceeds char vector size: {}",
                         handle.initialDataOffset(), mem_count_.charCount());
    return false;

  }

  return static_cast<bool>(char_memory_[handle.initialDataOffset()]);

}


std::vector<int64_t> kgl::DataMemoryBlock::getInteger(const InfoResourceHandle& handle) const {

  std::vector<int64_t> integer_vec;

  if (handle.dynamicType() == DataDynamicType::FixedData) {

    if (handle.initialDataOffset() >= mem_count_.integerCount()) {

      ExecEnv::log().error("DataMemoryBlock::getInteger, Handle offset: {} exceeds integer vector size: {}",
                           handle.initialDataOffset(), mem_count_.integerCount());
      return integer_vec;

    }

    InfoIntegerType scalar_integer = integer_memory_[handle.initialDataOffset()];

    if (scalar_integer != VCFInfoParser::MISSING_VALUE_INTEGER_) {

      integer_vec.push_back(scalar_integer);

    }

  } else if (handle.dynamicType() == DataDynamicType::FixedDynamic) {

    if (handle.initialDataOffset() >= mem_count_.integerCount()) {

      ExecEnv::log().error("DataMemoryBlock::getInteger, Handle offset: {} exceeds integer vector size: {}",
                           handle.initialDataOffset(), mem_count_.integerCount());
      return integer_vec;

    }

    InfoIntegerType scalar_integer = integer_memory_[handle.initialDataOffset()];

    // If a missing value then we assume it is stored as an array, however it could just be missing.
    if (scalar_integer == VCFInfoParser::MISSING_VALUE_INTEGER_) {

      std::optional<const InfoArrayIndex> array_index_opt = findArrayIndex(handle.handleId());

      if (array_index_opt) {

        const InfoArrayIndex array_index = array_index_opt.value();

        if (array_index.infoOffset() + array_index.infoSize() > mem_count_.integerCount()) {

          ExecEnv::log().error("DataMemoryBlock::getInteger, Handle: {},  FixedDynamic Array Offset: {} + Array Size: {} exceeds Integer Array Size: {}",
                               handle.handleId(), array_index.infoOffset(), array_index.infoSize(), mem_count_.integerCount());
          return integer_vec;

        }

        for (size_t index = 0; index < array_index.infoSize(); ++index) {

          integer_vec.push_back(integer_memory_[array_index.infoOffset() + index]);

        }

      }

    } else {  // Valid value so assume a scalar.

      integer_vec.push_back(scalar_integer);

    }

  } else if (handle.dynamicType() == DataDynamicType::DynamicData) {

    std::optional<const InfoArrayIndex> array_index_opt = findArrayIndex(handle.handleId());

    if (not array_index_opt) {

      ExecEnv::log().error("DataMemoryBlock::getInteger, Array index for data item handle: {} not found", handle.handleId());
      return integer_vec;

    }

    const InfoArrayIndex array_index = array_index_opt.value();

    if (array_index.infoOffset() + array_index.infoSize() > mem_count_.integerCount()) {

      ExecEnv::log().error("DataMemoryBlock::getInteger, Handle: {}, DynamicData Array Offset: {} + Array Size: {} exceeds Integer Array Size: {}",
                           handle.handleId(), array_index.infoOffset(), array_index.infoSize(), mem_count_.integerCount());
      return integer_vec;

    }

    for (size_t index = 0; index < array_index.infoSize(); ++index) {

      integer_vec.push_back(integer_memory_[array_index.infoOffset() + index]);

    }

  }

  return integer_vec;

}


std::vector<double> kgl::DataMemoryBlock::getFloat(const InfoResourceHandle& handle) const {

  std::vector<double> float_vec;

  if (handle.dynamicType() == DataDynamicType::FixedData) {

    if (handle.initialDataOffset() >= mem_count_.floatCount()) {

      ExecEnv::log().error("DataMemoryBlock::getFloat, Handle offset: {} exceeds float vector size: {}",
                           handle.initialDataOffset(), mem_count_.floatCount());
      return float_vec;

    }

    InfoFloatType scalar_float = float_memory_[handle.initialDataOffset()];

    if (scalar_float != VCFInfoParser::MISSING_VALUE_FLOAT_) {

      float_vec.push_back(scalar_float);

    }

  } else if (handle.dynamicType() == DataDynamicType::FixedDynamic) {

    if (handle.initialDataOffset() >= mem_count_.floatCount()) {

      ExecEnv::log().error("DataMemoryBlock::getFloat, Handle offset: {} exceeds float vector size: {}",
                           handle.initialDataOffset(), mem_count_.floatCount());
      return float_vec;

    }

    InfoFloatType scalar_float = float_memory_[handle.initialDataOffset()];

    // If a missing value then we assume it is stored as an array, however it could just be missing.
    if (scalar_float == VCFInfoParser::MISSING_VALUE_FLOAT_) {

      std::optional<const InfoArrayIndex> array_index_opt = findArrayIndex(handle.handleId());

      if (array_index_opt) {

        const InfoArrayIndex array_index = array_index_opt.value();

        if (array_index.infoOffset() + array_index.infoSize() > mem_count_.floatCount()) {

          ExecEnv::log().error("DataMemoryBlock::getFloat, Array Offset: {} + Array Size: {} exceeds Float Array Size: {}",
                               array_index.infoOffset(), array_index.infoSize(), mem_count_.floatCount());
          return float_vec;

        }

        for (size_t index = 0; index < array_index.infoSize(); ++index) {

          float_vec.push_back(float_memory_[array_index.infoOffset() + index]);

        }

      }

    } else {  // Valid value so assume a scalar.

      float_vec.push_back(scalar_float);

    }

  } else if (handle.dynamicType() == DataDynamicType::DynamicData) {

    std::optional<const InfoArrayIndex> array_index_opt = findArrayIndex(handle.handleId());

    if (not array_index_opt) {

      ExecEnv::log().error("DataMemoryBlock::getFloat, Array index for data item handle: {} not found", handle.handleId());
      return float_vec;

    }

    const InfoArrayIndex array_index = array_index_opt.value();

    if (array_index.infoOffset() + array_index.infoSize() > mem_count_.floatCount()) {

      ExecEnv::log().error("DataMemoryBlock::getFloat, Array Offset: {} + Array Size: {} exceeds Float Array Size: {}",
                           array_index.infoOffset(), array_index.infoSize(), mem_count_.floatCount());
      return float_vec;

    }

    for (size_t index = 0; index < array_index.infoSize(); ++index) {

      float_vec.push_back(float_memory_[array_index.infoOffset() + index]);

    }

  }

  return float_vec;

}

std::vector<std::string> kgl::DataMemoryBlock::getString(const InfoResourceHandle& handle) const {

  std::vector<std::string> string_vec;

  if (handle.dynamicType() == DataDynamicType::FixedData) {

    if (handle.initialDataOffset() >= mem_count_.stringCount()) {

      ExecEnv::log().error("DataMemoryBlock::getString, Handle offset: {} exceeds string vector size: {}",
                           handle.initialDataOffset(), mem_count_.stringCount());
      return string_vec;

    }

    std::string_view str_view = string_memory_[handle.initialDataOffset()];
    string_vec.emplace_back(std::string(str_view));

  } else if (handle.dynamicType() == DataDynamicType::DynamicData) {

    std::optional<const InfoArrayIndex> array_index_opt = findArrayIndex(handle.handleId());

    if (not array_index_opt) {

      ExecEnv::log().error("DataMemoryBlock::getString, Array index for data item handle: {} not found", handle.handleId());
      return string_vec;

    }

    const InfoArrayIndex array_index = array_index_opt.value();

    if (array_index.infoOffset() + array_index.infoSize() > mem_count_.stringCount()) {

      ExecEnv::log().error("DataMemoryBlock::getString, Array Offset: {} + Array Size: {} exceeds String Array Size: {}",
                           array_index.infoOffset(), array_index.infoSize(), mem_count_.stringCount());
      return string_vec;

    }

    for (size_t index = 0; index < array_index.infoSize(); ++index) {

      std::string_view str_view = string_memory_[array_index.infoOffset() + index];
      string_vec.emplace_back(std::string(str_view));

    }

  }

  return string_vec;

}

