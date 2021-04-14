//
// Created by kellerberrin on 28/02/18.
//

#ifndef KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H


#include "kel_utility.h"
#include "kgl_alphabet_dna5.h"
#include "kgl_alphabet_string.h"
#include "kgl_resource_db.h"
#include "kgl_variant_file_vcf_record.h"

#include <seqan/vcf_io.h>

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF seqan record parser.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseVCFRecord {

public:

  explicit ParseVCFRecord(const VcfRecord& vcf_record) : vcf_record_(vcf_record) {}
  ~ParseVCFRecord() = default;

  [[nodiscard]] bool parseRecord(const ContigId_t& contig_id, std::shared_ptr<const GenomeReference> genome_db_ptr);

  [[nodiscard]] size_t requiredFormatSize() const;
  [[nodiscard]] size_t GTOffset() const { return GT_offset_; }
  [[nodiscard]] size_t GQOffset() const { return GQ_offset_; }
  [[nodiscard]] size_t PLOffset() const { return PL_GL_offset_; }
  [[nodiscard]] size_t ADOffset() const { return AD_offset_; }
  [[nodiscard]] size_t DPOffset() const { return DP_offset_; }

  [[nodiscard]] Phred_t quality() const { return quality_; }

  [[nodiscard]] const std::string& reference() const { return reference_; }
  [[nodiscard]] const std::vector<std::string>& alleles() const { return alleles_; }

  [[nodiscard]] ContigOffset_t offset() const { return allele_offset_; }
  [[nodiscard]] std::shared_ptr<const ContigReference> contigPtr() const { return contig_opt_.value(); }

  [[nodiscard]] bool isSNP() const;

private:

  bool parse_result_;

  const VcfRecord& vcf_record_;
  std::vector<std::string> format_fields_;  // GT:AD:DP:GQ:PGT:PID:PL
  size_t GT_offset_;
  size_t GQ_offset_;
  size_t PL_GL_offset_;
  size_t AD_offset_;
  size_t DP_offset_;
  size_t required_size_;
  std::string reference_;
  std::vector<std::string> alleles_;
  ContigOffset_t allele_offset_;
  std::optional<std::shared_ptr<const ContigReference>> contig_opt_;
  Phred_t quality_;

  constexpr static const char* FORMAT_SEPARATOR_ = ":";
  constexpr static const char* ALLELE_SEPARATOR_ = ",";
  constexpr static const char* GT_ = "GT";
  constexpr static const char* GQ_ = "GQ";
  constexpr static const char* PL_ = "PL";    // Gatk uses PL (integers)
  constexpr static const char* GL_ = "GL";    // FreeBayes uses GL (floats)
  constexpr static const char* AD_ = "AD";
  constexpr static const char* DP_ = "DP";


  static void parseString(const std::string& parse_string,
                          const std::string& separator_string,
                          std::vector<std::string>& parse_results);

  static bool findString(const std::vector<std::string>& string_vec,
                         const std::string& search_string,
                         size_t& offset);


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF seqan genotype parser.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using FormatOffset = std::pair<size_t, size_t>;  // .first is the offset (in chars), .second is the size (in chars)
using FormatArray = std::array<FormatOffset, 100>; // The parsed formats

class ParseVCFGenotype {

public:

  ParseVCFGenotype(const std::string& format_char_string) : format_count_(0), parse_result_(false) {

    parseGenotype(format_char_string);

  }
  ~ParseVCFGenotype() = default;


  // Get a format field as a string.
  std::string getFormatString(size_t format_offset, const std::string& format_char_string) const;

  // Get the first char of a format field.
  char getFormatChar(size_t format_offset, const std::string& format_char_string) const;

  size_t formatCount() const { return format_count_; }

  const FormatArray& formatOffsets() const { return format_offsets_; }

  constexpr static const char FORMAT_SEPARATOR_ = ':';


private:

  constexpr static const size_t MAX_FORMAT_SIZE_ = 100;

  // Parse format fields
  bool parseGenotype(const std::string& format_char_string);

  size_t format_count_;
  FormatArray format_offsets_;
  bool parse_result_;


};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Diploid genotype.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DiploidAlleles = std::vector<std::pair<size_t, size_t>>;  // .first is A, .second is B.

class DiploidGenotypes {

public:

  DiploidGenotypes() = default;
  ~DiploidGenotypes() = default;

  std::string genotypeText(size_t allele_count) const;
  DiploidAlleles generateGenotype(size_t allele_count) const;
  void generateGenotypeVector(size_t max_alleles);

private:

  std::vector<DiploidAlleles> diploid_alleles_;




};


}   // end namespace



#endif //KGL_KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H
