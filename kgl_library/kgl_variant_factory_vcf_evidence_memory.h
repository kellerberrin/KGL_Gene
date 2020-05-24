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

  void unityArrayCount(size_t count) { unity_array_count_ = count; }
  void arrayCount(size_t count) { array_count_ = count; }
  void floatCount(size_t count) { float_count_ = count; }
  void integerCount(size_t count) { integer_count_ = count; }
  void stringCount(size_t count) { string_count_ = count; }
  void charCount(size_t count) { char_count_ = count; }

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
// Internal array index structure (only uses 64 bits for size efficiency).
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

  // Minimize storage size (64 bits).
  using VariableIndexImpl = uint16_t;
  using InfoOffsetImpl = uint32_t;
  using InfoCountImpl = uint16_t;

  VariableIndexImpl info_variable_index_{0};  // corresponds to the variable index in the header.
  InfoOffsetImpl info_element_offset_{0};
  InfoCountImpl info_element_count_{0};

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The data object contains only 5 data structures
// 1. A vector of Floats.
// 2. A vector of ints.
// 3. A vector of chars.
// 4  A vector of std::string_view to index into the char vector to define strings.
// 5. A vector of { offset, size, type } to define integer or float vectors at run time when parsed info is presented.
// 6. A vector of { offset, size, type } to define runtime vectors if the internal_unity_integer and internal_unity_float
// are not scalar.
// If we assume that AlternateAllele, AllAllele, Genotype, only contain scalars then 1 to 5 are known at subscribe time.
// These indexes are held for each subscribed variable below.
//
// If, however, the AlternateAllele, AllAllele, Genotype subscribers are vectors then 4. and 5. are calculated when the parsed INFO
// data is presented to the InfoDataBlock object and the data vectors for 1. to 5. are adjusted accordingly.

// To simplify initial coding, the AlternateAllele, AllAllele, Genotype are assumed to be scalars.
// If they are specified as array then this treated as infrequent case

// Forward declarations.
class InfoEvidenceHeader;
class InfoDataBlock {

public:

  InfoDataBlock(std::shared_ptr<const InfoEvidenceHeader> info_evidence_header)
  : info_evidence_header_(std::move(info_evidence_header)) {}
  InfoDataBlock(const InfoDataBlock &) = delete;
  ~InfoDataBlock() = default;

  [[nodiscard]] bool indexAndVerify(size_t field_address,     // The field index
                                    size_t field_id,           // the field identifier
                                    InfoEvidenceIntern internal_type,
                                    const std::optional<InfoParserToken> &token,
                                    InfoDataUsageCount &dynamic_accounting,
                                    InfoDataUsageCount &static_accounting);

  std::optional<InfoArrayIndex> findUnityArrayIndex(size_t field_index);

  void allocateMemory(const InfoDataUsageCount &type_count);
  [[nodiscard]] const InfoDataUsageCount &getRawMemoryUsage() const { return type_count_; }

  // Constants for missing values.
  // This is a bit dubious because these are potentially (though unlikely) valid field values.
  // However the alternative is to setup a system of missing value flags for these field types.
  constexpr static const InfoParserFloat MISSING_VALUE_FLOAT = std::numeric_limits<InfoParserFloat>::lowest();
  constexpr static const InfoParserInteger MISSING_VALUE_INTEGER = std::numeric_limits<InfoParserInteger>::lowest();


private:


  std::shared_ptr<const InfoEvidenceHeader> info_evidence_header_; // The data header and indexing structure.
  InfoDataUsageCount type_count_;  // Size count object.

  // Raw memory to efficiently store the info data.
  std::unique_ptr<char[]> char_memory_;
  std::unique_ptr<InfoIntegerType[]> integer_memory_;
  std::unique_ptr<InfoFloatType[]> float_memory_;
  std::unique_ptr<InfoArrayIndex[]> array_memory_;
  std::unique_ptr<InfoArrayIndex[]> unity_array_memory_;
  std::unique_ptr<std::string_view[]> string_memory_;


  // Implementation functions of indexAndVerify for each internal data type.
  // These exist to breakup blocks of unsightly code into small manageable chunks.
  bool internChar(size_t field_address,
                  const std::optional<InfoParserToken> &token,
                  InfoDataUsageCount &static_accounting);

  bool internInteger(size_t field_address,
                      const std::optional<InfoParserToken> &token,
                     InfoDataUsageCount &static_accounting);

  bool internFloat(size_t field_address,
                   const std::optional<InfoParserToken> &token,
                   InfoDataUsageCount &static_accounting);

  bool internString(size_t field_address,
                    const std::optional<InfoParserToken> &token,
                    InfoDataUsageCount &dynamic_accounting,
                    InfoDataUsageCount &static_accounting);

  bool internIntegerArray(size_t field_address,
                          size_t field_id,
                          const std::optional<InfoParserToken> &token,
                          InfoDataUsageCount &dynamic_accounting,
                          InfoDataUsageCount &static_accounting);

  bool internFloatArray(size_t field_address,
                        size_t field_id,
                        const std::optional<InfoParserToken> &token,
                        InfoDataUsageCount &dynamic_accounting,
                        InfoDataUsageCount &static_accounting);

  bool internStringArray(size_t field_address,
                         size_t field_id,
                         const std::optional<InfoParserToken> &token,
                         InfoDataUsageCount &dynamic_accounting,
                         InfoDataUsageCount &static_accounting);

  bool internUnityInteger(size_t field_offset,
                          size_t field_index,
                          const std::optional<InfoParserToken> &token,
                          InfoDataUsageCount &dynamic_accounting,
                          InfoDataUsageCount &static_accounting);

  bool internUnityFloat(size_t field_offset,
                        size_t field_index,
                        const std::optional<InfoParserToken> &token,
                        InfoDataUsageCount &dynamic_accounting,
                        InfoDataUsageCount &static_accounting);

  // Implementation functions that retrieve data for each internal data type.

  std::optional<bool> getBoolean(size_t field_address);

  std::optional<InfoIntegerType> internInteger(size_t field_address);

  std::optional<InfoFloatType> internFloat(size_t field_address);

  std::optional<std::string_view> internString(size_t field_address);

  InfoArrayIndex internIntegerArray(size_t field_address, size_t field_id);

  InfoArrayIndex internFloatArray(size_t field_address, size_t field_id);

  InfoArrayIndex internStringArray(size_t field_address, size_t field_id);

  InfoArrayIndex internUnityInteger(size_t field_offset, size_t field_index);

  InfoArrayIndex internUnityFloat(size_t field_offset, size_t field_index);

};



} //namespace.


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEMORY_H
