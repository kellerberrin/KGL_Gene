//
// Created by kellerberrin on 24/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_BLOCK_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_BLOCK_H



#include "kel_exec_env.h"
#include "kgl_variant_factory_vcf_evidence_data.h"
#include "kgl_variant_factory_vcf_evidence_memory.h"


namespace kellerberrin::genome {   //  organization level namespace

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

  void allocateMemory(const InfoDataUsageCount &type_count);
  [[nodiscard]] const InfoDataUsageCount &getRawMemoryUsage() const { return type_count_; }

  // Constants for missing values.
  // This is a bit dubious because these are potentially (though unlikely) valid field values.
  // However the alternative is to setup a system of missing value flags for these field types.
  constexpr static const InfoParserFloat MISSING_VALUE_FLOAT_ = std::numeric_limits<InfoParserFloat>::lowest();
  constexpr static const InfoParserInteger MISSING_VALUE_INTEGER_ = std::numeric_limits<InfoParserInteger>::lowest();

  // Implementation functions that retrieve data for each internal data type.

  size_t getDataSize(size_t field_address,     // The field index
                     size_t field_id,           // the field identifier
                     InfoEvidenceIntern internal_type) const;

  bool getBoolean(size_t field_address) const;
  InfoIntegerType getInteger(size_t field_address) const;
  InfoFloatType getFloat(size_t field_address) const;
  std::string getString(size_t field_address) const;
  std::vector<int64_t> getIntegerArray(size_t field_address, size_t field_id) const;
  std::vector<double> getFloatArray(size_t field_address, size_t field_id) const;
  std::vector<std::string> getStringArray(size_t field_address, size_t field_id) const;
  std::vector<int64_t> getUnityIntegerArray(size_t field_offset, size_t field_id) const;
  std::vector<double> getUnityFloatArray(size_t field_offset, size_t field_id) const;

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

// Data accessor functions. No data is copied used these functions.
  std::string_view getStringView(size_t field_address) const;
  InfoArrayIndex getArray(size_t field_address, size_t field_id) const;
  InfoArrayIndex getIntegerUnityArray(size_t field_address, size_t field_id) const;
  InfoArrayIndex getFloatUnityArray(size_t field_address, size_t field_id) const;

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


};



} // name space

#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_BLOCK_H
