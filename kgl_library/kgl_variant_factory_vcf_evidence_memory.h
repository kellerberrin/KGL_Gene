//
// Created by kellerberrin on 23/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H

#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_vcf_evidence_data.h"



namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object manages the storage of individual Info fields from the raw data allocated below.

struct InfoDataUsageCount {

public:

  InfoDataUsageCount() = default;
  ~InfoDataUsageCount() = default;

  bool operator==(const InfoDataUsageCount& cmp) const;

  [[nodiscard]] size_t unityArrayCount() const { return unity_array_count_; }
  [[nodiscard]] size_t arrayCount() const { return array_count_; }
  [[nodiscard]] size_t floatCount() const { return float_count_; }
  [[nodiscard]] size_t integerCount() const { return integer_count_; }
  [[nodiscard]] size_t stringCount() const { return string_count_; }
  [[nodiscard]] size_t charCount() const { return char_count_; }

  void unityArrayCountAdd(size_t count) { unity_array_count_ += count; }
  void arrayCountAdd(size_t count) { array_count_ += count; }
  void floatCountAdd(size_t count) { float_count_ += count; }
  void integerCountAdd(size_t count) { integer_count_ += count; }
  void stringCountAdd(size_t count) { string_count_ += count; }
  void charCountAdd(size_t count) { char_count_ += count; }

  // Notionally allocates data on the 5 data arrays and returns a data block
  // The function is to be used sequentially on all subscribed variables.
  // The final values are used to actually allocate memory in the InfoDataBlock object.

  // Set up the indexes and pre-allocate fixed sized fields (run once for all fields)
  [[nodiscard]] size_t staticIncrementAndAllocate(InfoEvidenceIntern internal_type);

  // Allocate additional memory space at runtime run for every Info field parsed.
  // Sets up all the array indexes and verifies the size of the total allocated memory.
  [[nodiscard]] bool dynamicIncrementAndAllocate(InfoEvidenceIntern internal_type, const InfoParserToken &token);

private:

  size_t unity_array_count_{0};
  size_t array_count_{0};
  size_t float_count_{0};
  size_t integer_count_{0};
  size_t string_count_{0};
  size_t char_count_{0};

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Internal array index structure (can use 64 bits for size efficiency).

// An array block, 1 for each data array. VCF data arrays are generally small, so this can be a memory overhead.
// Small is 64 bits total, large is 2 x 64 = 128 bits, generally large is fine (unless space is tight and a lot of arrays).
//#define INFO_DATA_SMALL_ARRAY_SIZE 1  // Uncomment for 64 bits per array, else 128 bits.


struct InfoArrayIndex {

public:

  InfoArrayIndex() = default;
  InfoArrayIndex(size_t variable_index, size_t element_offset, size_t element_count) {

    infoVariableIndex(variable_index);
    infoOffset(element_offset);
    infoSize(element_count);

  }
  InfoArrayIndex(const InfoArrayIndex&) = default;
  ~InfoArrayIndex() = default;

  [[nodiscard]] size_t infoVariableIndex() const { return static_cast<size_t>(info_variable_index_); }
  [[nodiscard]] size_t infoOffset() const { return static_cast<size_t>(info_element_offset_); }
  [[nodiscard]] size_t infoSize() const { return static_cast<size_t>(info_element_count_); }

  void infoVariableIndex(size_t info_variable_index) { info_variable_index_ = static_cast<VariableIndexImpl>(info_variable_index); }
  void infoOffset(size_t info_element_offset) { info_element_offset_ = static_cast<InfoOffsetImpl>(info_element_offset); }
  void infoSize(size_t info_element_count) { info_element_count_ = static_cast<InfoCountImpl>(info_element_count); }


private:

#ifdef INFO_DATA_SMALL_ARRAY_SIZE
  // Minimize storage size (64 bits).
  using VariableIndexImpl = uint16_t;
  using InfoOffsetImpl = uint32_t;
  using InfoCountImpl = uint16_t;
#else
  using VariableIndexImpl = uint32_t;
  using InfoOffsetImpl = uint64_t;
  using InfoCountImpl = uint32_t;
#endif


  VariableIndexImpl info_variable_index_{0};  // corresponds to the variable index in the header.
  InfoOffsetImpl info_element_offset_{0};
  InfoCountImpl info_element_count_{0};

};



} //namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H
