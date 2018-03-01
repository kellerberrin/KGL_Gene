//
// Created by kellerberrin on 28/02/18.
//

#include "kgl_variant_factory_record_vcf_impl.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF seqan record parser.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ParseVCFRecord::parseRecord(const ContigId_t& contig_id, std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  bool parse_result = true;

  // Get the format fields for Genetype analysis.
  parseString(seqan::toCString(vcf_record_.format), FORMAT_SEPARATOR_, format_fields_);

  // Get the GT format offset.
  if (not findString(format_fields_, GT_, GT_offset_)) {

    ExecEnv::log().error("Format field: {} not found", GT_);
    parse_result = false;

  }

  // Get the GQ format offset.
  if (not findString(format_fields_, GQ_, GQ_offset_)) {

    ExecEnv::log().error("Format field: {} not found", GQ_);
    parse_result = false;

  }

  // Get the PL format offset.
  if (not findString(format_fields_, PL_, PL_offset_)) {

    ExecEnv::log().error("Format field: {} not found", PL_);
    parse_result = false;

  }

  required_size_ = requiredFormatSize();  // Must have this many format fields in a Genotype.

  // Get the reference DNA
  reference_ = seqan::toCString(vcf_record_.ref);

  // Get the allelle DNA vector.
  parseString(seqan::toCString(vcf_record_.alt), ALLELE_SEPARATOR_, alleles_);

  // Get the offset.
  allele_offset_ = static_cast<ContigOffset_t >(vcf_record_.beginPos);

  // Get the contig pointer.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().error("Contig: {} is not in the Genome Database", contig_id);
    contig_ptr = nullptr;
    parse_result = false;

  }

  std::shared_ptr<const DNA5SequenceLinear> contig_ref = contig_ptr->sequence().unstrandedRegion(allele_offset_, reference_.length());

  if (contig_ref->getSequenceAsString() != reference_) {

    ExecEnv::log().error("Variant reference: {} does not match Contig region: {} at offset: {}",
                         reference_, contig_ref->getSequenceAsString(), allele_offset_);
    parse_result = false;

  }

  // Get the overall allelle quality
  quality_ = static_cast<Phred_t>(vcf_record_.qual);

  parse_result_ = parse_result;

  return parse_result;

}


void kgl::ParseVCFRecord::parseString(const std::string& parse_string,
                                      const std::string& separator_string,
                                      std::vector<std::string>& parse_results) {

  bt::char_separator<char> item_key_sep(separator_string.c_str());
  bt::tokenizer<bt::char_separator<char>> tokenize_item(parse_string, item_key_sep);
  for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

    parse_results.push_back(*iter_item);

  }

}


bool kgl::ParseVCFRecord::findString(const std::vector<std::string>& string_vec,
                                     const std::string& search_string,
                                     size_t& offset) {

  bool found_string = false;
  offset = 0;

  for (const auto& string_item : string_vec) {

    if (string_item == search_string) {

      found_string = true;
      break;

    }

    ++offset;

  }

  return found_string;

}


size_t kgl::ParseVCFRecord::requiredFormatSize() const {

  size_t required_size = 0;

  required_size = GT_offset_ >= GQ_offset_ ? GT_offset_ : GQ_offset_;
  required_size = PL_offset_ >= required_size ? PL_offset_ : required_size;
  ++required_size; // Convert offset to a size.

  return required_size;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF seqan genotype parser.
// .first is the offset (in chars), .second is the size (in chars)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ParseVCFGenotype::parseGenotype(const seqan::CharString& format_char_string) {

  bool parse_result = true;
  size_t record_index = 0;
  size_t field_index = 0;
  size_t field_size = 0;


  for (size_t idx = 0; idx < seqan::length(format_char_string); ++idx)
  {

    if (format_char_string[idx] == FORMAT_SEPARATOR_) {

      format_offsets_[record_index].first = field_index; // index
      format_offsets_[record_index].second = field_size; // size;
      ++record_index;
      field_index = idx;
      ++field_index;
      field_size = 0;

    } else {

      ++field_size;

    }

  }

  if (field_size > 0) {

    format_offsets_[record_index].first = field_index; // index
    format_offsets_[record_index].second = field_size; // size;
    format_count_ = record_index;
    ++format_count_;

  } else {

    format_count_ = record_index;

  }

  parse_result_ = parse_result;

  return parse_result;

}


std::string kgl::ParseVCFGenotype::getFormatString(size_t format_offset, const seqan::CharString &format_char_string) const {

  if (!parse_result_) {

    ExecEnv::log().error("ParseVCFGenotype(), Genotype has not been parsed");
    return "";

  }

  if (format_offset >= format_count_) {

    ExecEnv::log().error("getFormatChar(), format offset: {} out of range, format count: {}", format_offset, format_count_);
    return "";

  }

  size_t ptr_offset = formatOffsets()[format_offset].first;
  size_t size = formatOffsets()[format_offset].second;
  const char* char_ptr = &(seqan::toCString(format_char_string)[ptr_offset]);

  std::string format_string(char_ptr, size);

  return format_string;

}


char kgl::ParseVCFGenotype::getFormatChar(size_t format_offset, const seqan::CharString &format_char_string) const {

  if (!parse_result_) {

    ExecEnv::log().error("ParseVCFGenotype(), Genotype has not been parsed");
    return ' ';

  }

  if (format_offset >= format_count_) {

    ExecEnv::log().error("getFormatChar(), format offset: {} out of range, format count: {}", format_offset, format_count_);
    return ' ';

  }

  size_t first_PL_char_offset = formatOffsets()[format_offset].first;

  return format_char_string[first_PL_char_offset];

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF Diploid genotypes for the format PL field.
// .first = A, .second = B
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DiploidAlleles kgl::DiploidGenotypes::generateGenotype(size_t allele_count) {

  DiploidAlleles diploid_alleles;

  for (size_t j = 0; j <= allele_count; ++j) {

    for (size_t k = 0; k <= j; ++k) {

      diploid_alleles.push_back(std::pair<size_t, size_t>(k, j));

    }

  }

  return diploid_alleles;

}


void kgl::DiploidGenotypes::generateGenotypeVector(size_t max_alleles) {

  for(size_t allele_count = 1; allele_count <= max_alleles; ++allele_count) {

    diploid_alleles_.push_back(generateGenotype(allele_count));

  }

}


std::string kgl::DiploidGenotypes::genotypeText(size_t allele_count) const {


  if (allele_count > diploid_alleles_.size() or allele_count < 1) {

    ExecEnv::log().error("genotypeText(), allele count: {} out of range, max alleles defined: {}",
                         allele_count , diploid_alleles_.size());
    return "";

  }

  --allele_count; // convert to a vector index offset.

  const DiploidAlleles& diploid_allele = diploid_alleles_[allele_count];

  std::stringstream allele_ss;
  for (auto allele_item : diploid_allele) {

    allele_ss << "(A:" << allele_item.first << ", B:" << allele_item.second << "), ";

  }

  return allele_ss.str();

}

