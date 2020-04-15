//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H


#include <string>
#include <vector>

#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Parses the VCF header for Contig and Genome infromation.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  bool parseHeader(const std::string& file_name);

  [[nodiscard]] const std::vector<std::string>& getGenomes() const { return vcf_genomes_; }
  [[nodiscard]] const std::vector<std::pair<ContigId_t, ContigSize_t>>& getContigs() const { return contig_info_; }

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record
  std::vector<std::pair<ContigId_t, ContigSize_t>> contig_info_;  // Contigs and lengths in the VCF file.

  // Parser constants.
  static constexpr const char* CONTIG_NAME_FRAGMENT_{"##contig"};
  static constexpr const size_t CONTIG_NAME_FRAGMENT_LENGTH_{8};
  static constexpr const char* CONTIG_INFO_START_{"=<"};
  static constexpr const char* CONTIG_INFO_END_{">"};
  static constexpr const char* FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.


};


} // end namespace


#endif //KGL_KGL_VARIANT_FILE_H
