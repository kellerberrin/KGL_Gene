//
// Created by kellerberrin on 28/02/18.
//

#ifndef KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H


#include "kgl_utility.h"
#include "kgl_alphabet_dna5.h"
#include "kgl_alphabet_string.h"
#include "kgl_genome_db.h"

#include <seqan/vcf_io.h>

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF seqan record parser.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseVCFRecord {

public:

  ParseVCFRecord(const seqan::VcfRecord& vcf_record,
                 const ContigId_t& contig_id,
                 std::shared_ptr<const GenomeDatabase> genome_db_ptr) : vcf_record_(vcf_record) {

    parse_result_ = parseRecord(contig_id, genome_db_ptr);

  }
  ~ParseVCFRecord() = default;

  size_t requiredFormatSize() const;
  size_t GTOffset() const { return GT_offset_; }
  size_t GQOffset() const { return GQ_offset_; }
  size_t PLOffset() const { return PL_offset_; }

private:

  bool parse_result_;

  const seqan::VcfRecord& vcf_record_;
  std::vector<std::string> format_fields_;  // GT:AD:DP:GQ:PGT:PID:PL
  size_t GT_offset_;
  size_t GQ_offset_;
  size_t PL_offset_;
  size_t required_size_;
  AlphabetString<DNA5> reference_;
  std::vector<std::string> alleles_;
  ContigOffset_t allele_offset_;
  std::shared_ptr<const ContigFeatures> contig_ptr_;
  Phred_t quality_;

  constexpr static const char* FORMAT_SEPARATOR_ = ":";
  constexpr static const char* ALLELE_SEPARATOR_ = ",";
  constexpr static const char* GT_ = "GT";
  constexpr static const char* GQ_ = "GQ";
  constexpr static const char* PL_ = "PL";

  bool parseRecord(const ContigId_t& contig_id, std::shared_ptr<const GenomeDatabase> genome_db_ptr);

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

  ParseVCFGenotype() : format_count_(0), parse_result_(false) {}
  ~ParseVCFGenotype() = default;

  bool parseGenotype(const seqan::CharString& format_char_string);

  std::string getPLstring(size_t PLoffset, const seqan::CharString& format_char_string) const;

  size_t formatCount() const { return format_count_; }

  const FormatArray& formatOffsets() const { return format_offsets_; }

private:


  constexpr static const size_t MAX_FORMAT_SIZE_ = 100;
  constexpr static const char FORMAT_SEPARATOR_ = ':';

  size_t format_count_;
  FormatArray format_offsets_;
  bool parse_result_;


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H
