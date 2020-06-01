//
// Created by kellerberrin on 30/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_BLK_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_BLK_H

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
// If we assume that AlternateAllele, AllAllele, Genotype, only contain scalars then 1 to 5 are known at subscribe time.
// These indexes are held for each subscribed variable below.
//
// If, however, the AlternateAllele, AllAllele, Genotype subscribers are vectors then 4. and 5. are calculated when the parsed INFO
// data is presented to the InfoDataMem object and the data vectors for 1. to 5. are adjusted accordingly.

// To minimize data allocation, the AlternateAllele, AllAllele, Genotype are assumed to be scalars.
// If they are specified as array then this is treated as an infrequent case

// Forward declarations.
class InfoEvidenceHeader;
class DataMemoryBlock {

public:

  DataMemoryBlock( std::shared_ptr<const InfoEvidenceHeader> info_evidence_header,
                   const InfoMemoryResource& initial_memory_resource,
                   const VCFInfoParser& info_parser);
  DataMemoryBlock(const DataMemoryBlock &) = delete;
  ~DataMemoryBlock() = default;

  // Constants for missing values.
  // This is a bit dubious because these are potentially (though unlikely) valid field values.
  constexpr static const InfoParserFloat MISSING_VALUE_FLOAT_ = std::numeric_limits<InfoParserFloat>::lowest();
  constexpr static const InfoParserInteger MISSING_VALUE_INTEGER_ = std::numeric_limits<InfoParserInteger>::lowest();

  // Implementation functions that retrieve data for each internal data type.
  bool getBoolean(const InfoResourceHandle& handle) const;
  std::vector<int64_t> getInteger(const InfoResourceHandle& handle) const;
  std::vector<double> getFloat(const InfoResourceHandle& handle) const;
  std::vector<std::string> getString(const InfoResourceHandle& handle) const;

  [[nodiscard]] const MemDataUsage& getUsageCount() const { return mem_count_; }

private:


  std::shared_ptr<const InfoEvidenceHeader> info_evidence_header_; // The data header and indexing structure.
  MemDataUsage mem_count_;  // Size count object.

  // Raw memory to efficiently store the info data.
  std::unique_ptr<char[]> char_memory_;
  std::unique_ptr<InfoIntegerType[]> integer_memory_;
  std::unique_ptr<InfoFloatType[]> float_memory_;
  std::unique_ptr<InfoArrayIndex[]> array_memory_;
  std::unique_ptr<std::string_view[]> string_memory_;

  // Lookup the array index.
  [[nodiscard]] std::optional<const InfoArrayIndex> findArrayIndex(size_t identifier) const;

  // Data write functions. Write the content of a parser token into the data block.
  void storeBoolean(const InfoMemoryResource& memory_resource, const InfoResourceHandle& handle, std::optional<const InfoParserToken> token);
  void storeString(const InfoResourceHandle& handle, std::optional<const InfoParserToken> token, FixedResourceInstance& string_usage);
  void storeInteger(const InfoResourceHandle& handle, std::optional<const InfoParserToken> token);
  void storeFloat(const InfoResourceHandle& handle, std::optional<const InfoParserToken> token);

};



} // name space


#endif //KGL_KGL_VARIANT_FACTORY_VCF_EVIDENCE_DATA_BLK_H
