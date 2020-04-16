//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H


#include <string>
#include <vector>

#include "kgl_genome_types.h"

#include <seqan/vcf_io.h>


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



class VcfRecord
{
public:

  // Default constructor.
  VcfRecord() : offset(INVALID_POS), qual(MISSING_QUAL) {}
  VcfRecord(seqan::VcfRecord&& vcf_record, ContigId_t&& contig_id);

  // Numeric id of the reference sequence.
  ContigId_t contig_id;
  // Position on the reference.
  ContigOffset_t offset;
  // Textual identifier of the variant.
  std::string id;
  // Bases in the reference.
  std::string ref;
  // Bases in the alternatives, comma-separated.
  std::string alt;
  // Quality
  double qual;
  // Value of FILTER field.
  std::string filter;
  // Value of INFO field.
  std::string info;
  // Value of FORMAT field.
  std::string format;
  // The genotype infos.
  std::vector<std::string> genotypeInfos;

private:
  // Constant for invalid position.
  static constexpr const ContigOffset_t INVALID_POS = std::numeric_limits<ContigOffset_t>::max();
  // Undefined quality number.
  static constexpr const double MISSING_QUAL = std::numeric_limits<double>::lowest();

};



} // end namespace


#endif //KGL_KGL_VARIANT_FILE_H
