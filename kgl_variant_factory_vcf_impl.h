//
// Created by kellerberrin on 29/12/17.
//

#ifndef KGL_VARIANT_FACTORY_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_VCF_IMPL_H



#include "kgl_utility.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_single.h"
#include "kgl_variant_factory_single.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <seqan/vcf_io.h>

namespace kgl = kellerberrin::genome;
namespace bt = boost;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory::VcfFileImpl does all the heavy lifting using 3rd a party library. In this case; Seqan.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ActiveContigMap = std::map<kgl::ContigId_t, kgl::ContigSize_t>;

class kgl::VcfFactory::VcfFileImpl {

public:

  VcfFileImpl() = default;
  ~VcfFileImpl() = default;

  std::shared_ptr<GenomeVariant> readParseVcfFile(const std::string& genome_name,
                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                  const std::string& vcf_file_name,
                                                  Phred_t variant_quality);

private:

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";
  constexpr static const char* ID_READ_DEPTH_ = "DPB";
  constexpr static const char* ID_PROPORTION_ = "AF";

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 5000;

  size_t vcf_record_count_;
  size_t vcf_record_ignored_;
  size_t vcf_record_error_;
  size_t vcf_record_rejected_;

  bool parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                      const seqan::VcfHeader& header,
                      ActiveContigMap& active_contig_map) const;

  bool parseVcfRecord(const std::string& genome_name,
                      const seqan::VcfRecord& record,
                      ActiveContigMap& active_contig_map,
                      std::shared_ptr<const ContigFeatures> contig_ptr,
                      std::shared_ptr<GenomeVariant> genome_variants,
                      Phred_t variant_quality,
                      bool& quality_ok) const;

  // Parse 1M ... XM in the cigar.
  bool parseCheck(size_t cigar_count,
                  std::shared_ptr<const ContigFeatures> contig_ptr,
                  const std::string& reference,
                  const std::string& alternate,
                  size_t& reference_index,
                  size_t& alternate_index,
                  ContigOffset_t& contig_offset) const;

  // Parse 1X ... XX in the cigar.

  bool parseSNP(size_t cigar_count,
                const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                std::shared_ptr<GenomeVariant> genome_variants,
                Phred_t quality,
                const std::string& info,
                const std::string& reference,
                const std::string& alternate,
                size_t& reference_index,
                size_t& alternate_index,
                ContigOffset_t& contig_offset,
                size_t& variant_count) const;

  // Parse 1I ... XI in the cigar.
  bool parseInsert(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& alternate,
                   size_t& alternate_index,
                   ContigOffset_t& contig_offset,
                   size_t& variant_count) const;

  // Parse 1D ... XD in the cigar.
  bool parseDelete(size_t cigar_count,
                   const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   size_t& reference_index,
                   ContigOffset_t& contig_offset,
                   size_t& variant_count) const;

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  bool tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                  std::map<std::string, std::string>& key_value_pairs) const;

// assumes input "key_1=value_1; ...;key_n=value_n"
  bool tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                std::map<std::string, std::string>& key_value_map) const;

  bool parseCigar(const std::string& cigar,
                  size_t& check_reference_size,
                  size_t& check_alternate_size,
                  std::vector<std::pair<char, size_t>>& parsed_cigar) const;

};



#endif //READSAMFILE_KGL_VARIANT_FACTORY_VCF_IMPL_H
