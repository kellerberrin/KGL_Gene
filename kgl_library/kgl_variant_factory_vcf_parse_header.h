//
// Created by kellerberrin on 11/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_PARSE_HEADER_H
#define KGL_VARIANT_FACTORY_VCF_PARSE_HEADER_H

#include "kgl_genome_types.h"
#include "kgl_genome_db.h"
#include "kgl_variant_file_impl.h"

#include <memory>
#include <string>
#include <vector>


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Structure to return VCF header information as a vector of pair<key, value>
using VcfHeaderInfo = std::vector<std::pair<std::string, std::string>>;
// A map of contigs.
using ActiveContigMap = std::map<ContigId_t, ContigSize_t>;
// A map of VCF INFO fields.
using VcfInfoKeyMap = std::map<std::string, std::string>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parses the VCF header for Contig and Genome information.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  [[nodiscard]] bool parseHeader(const std::string &file_name);

  [[nodiscard]] const std::vector<std::string> &getGenomes() const { return vcf_genomes_; }

  [[nodiscard]] const VcfHeaderInfo &getHeaderInfo() const { return vcf_header_info_; }

  [[nodiscard]] static bool parseVcfHeader( std::shared_ptr<const GenomeReference> genome_db_ptr,
                                            const VcfHeaderInfo& header,
                                            ActiveContigMap& active_contig_map,
                                            bool cigar_required);

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record
  VcfHeaderInfo vcf_header_info_;
  // Parser constants.
  static constexpr const char *KEY_SEPARATOR_{"="};
  static constexpr const char *KEY_PREFIX_{"##"};
  static constexpr const char *FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* ID_CIGAR_VALUE_ = "CIGAR";

  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  [[nodiscard]] static bool tokenizeVcfHeaderKeyValues( const std::string& key_value_text,
                                                        VcfInfoKeyMap& key_value_pairs);

};


} // namespace.


#endif //KGL_VARIANT_FACTORY_VCF_PARSE_HEADER_H
