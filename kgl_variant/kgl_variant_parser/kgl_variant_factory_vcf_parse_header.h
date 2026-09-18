//
// Created by kellerberrin on 11/5/20.
//

#ifndef KGL_VARIANT_FACTORY_VCF_PARSE_HEADER_H
#define KGL_VARIANT_FACTORY_VCF_PARSE_HEADER_H

#include "kgl_genome_types.h"
#include "kgl_genome_genome.h"

#include <memory>
#include <string>
#include <vector>


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Structure to return VCF header information as a vector of pair<key, value>
using VCFHeaderInfo = std::vector<std::pair<std::string, std::string>>;
// A map of contigs.
using VCFContigMap = std::map<ContigId_t, ContigSize_t>;
// A map of VCF INFO fields.
using VcfInfoKeyMap = std::map<std::string, std::string>;
// A remapping of VCF contig names to Reference Genome contig_ref_ptr names.
// The key is the VCF contig_ref_ptr ID, value is the Reference Genome ID.
// Matched by contig_ref_ptr sizes so there is a (small) chance of a mis-mapping
using VCFContigAliasMap = std::map<std::string, std::string>;

// A VCF INFO record.
struct VCFInfoRecord {

  std::string ID;
  std::string description;
  std::string type;
  std::string number;
  std::string source;
  std::string version;

};
// A VCF INFO record map.
using VCFInfoRecordMap = std::map<std::string, VCFInfoRecord>;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parses the VCF header for Contig and Genome information.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  [[nodiscard]] bool parseHeader(const std::string &file_name);

  [[nodiscard]] const std::vector<std::string> &getGenomes() const { return vcf_genomes_; }

  [[nodiscard]] const VCFHeaderInfo &getHeaderInfo() const { return vcf_header_info_; }

  [[nodiscard]] static bool parseVcfHeader(const VCFHeaderInfo& header, VCFContigMap& vcf_contig_map, VCFInfoRecordMap& vcf_info_map);

  [[nodiscard]] static bool checkVCFReferenceContigs(const VCFContigMap& vcf_contig_map, std::shared_ptr<const GenomeReference> reference_genome);

  [[nodiscard]] static bool VCFContigAliasRemapping(const VCFContigMap& vcf_contig_map,
                                                    std::shared_ptr<const GenomeReference> reference_genome,
                                                    VCFContigAliasMap contig_alias_map);

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record
  VCFHeaderInfo vcf_header_info_;

  // Parser constants.
  static constexpr const char *KEY_SEPARATOR_{"="};
  static constexpr const char *KEY_PREFIX_{"##"};
  static constexpr const char *FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const char RECORD_FIELD_LIST_SEPARATOR_{'\t'};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.

  constexpr static const char* HEADER_CONTIG_KEY_ = "CONTIG";
  constexpr static const char* ID_KEY_ = "ID";
  constexpr static const char* CONTIG_LENGTH_KEY_ = "LENGTH";
  constexpr static const char* HEADER_INFO_KEY_ = "INFO";
  constexpr static const char* DESCRIPTION_KEY_ = "DESCRIPTION";
  constexpr static const char* TYPE_KEY_ = "TYPE";
  constexpr static const char* NUMBER_KEY_ = "NUMBER";
  constexpr static const char* SOURCE_KEY_ = "SOURCE";
  constexpr static const char* VERSION_KEY_ = "VERSION";


  // assumes input "<key_1=value_1, ...,key_n=value_n>"
  [[nodiscard]] static bool tokenizeVcfHeaderKeyValues( const std::string& key_value_text,
                                                        VcfInfoKeyMap& key_value_pairs);

};


} // namespace.


#endif //KGL_VARIANT_FACTORY_VCF_PARSE_HEADER_H
