//
// Created by kellerberrin on 28/02/18.
//

#include "kgl_variant_factory_record_vcf_impl.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF seqan record parser.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::ParseVCFRecord::ParseVCFRecord( const std::string& genome_contig,
                                     const VCFRecord& vcf_record,
                                     const std::shared_ptr<const GenomeReference>& genome_db_ptr) {

  // Get the format fields for Genetype analysis.
  format_fields_ = Utility::charTokenizer(vcf_record.format, FORMAT_SEPARATOR_);

  // Get the reference DNA
  reference_ = vcf_record.ref;

  // Get the allelle DNA vector.
  alleles_ = Utility::charTokenizer(vcf_record.alt, ALLELE_SEPARATOR_);

  // Get the offset.
  allele_offset_ = vcf_record.offset;

  // Get the contig_ref_ptr pointer.
  auto contig_opt = genome_db_ptr->getContigSequence(genome_contig);
  if (not contig_opt) {

    ExecEnv::log().error("ParseVCFRecord::parseRecord; Contig: {} is not in the Genome Database", vcf_record.contig_id);
    parse_result_ =false;
    return;

  }

  contig_ptr_ = contig_opt.value();

  OpenRightUnsigned referenceInterval(allele_offset_, allele_offset_+reference_.length());
  auto contig_ref_opt = contig_ptr_->sequence().subSequence(referenceInterval);
  if (not contig_ref_opt) {

    ExecEnv::log().error("Cannot extract reference interval: {} from contig: {} contig_ref_ptr interval: {}",
                         referenceInterval.toString(),
                         contig_ptr_->contigId(),
                         contig_ptr_->sequence().interval().toString());
    parse_result_ = false;
    return;

  }
  const DNA5SequenceLinear& contig_ref = contig_ref_opt.value();
  if (contig_ref.getSequenceAsString() != reference_) {

    ExecEnv::log().error("Variant reference: {} does not match Contig region: {} at offset: {}",
                         reference_, contig_ref.getSequenceAsString(), allele_offset_);
    parse_result_ = false;
    return;

  }

  // Get the overall allele quality
  quality_ = static_cast<Phred_t>(vcf_record.qual);

  // Look at the filter field for "Pass"
  passed_filter_ = Utility::toupper(vcf_record.filter) == PASSED_FILTERS_;

}



std::optional<size_t> kgl::ParseVCFRecord::formatIndex(const std::string& format_code) const {

  bool found_string{false};
  size_t offset{0};

  for (const auto& string_item : format_fields_) {

    if (string_item == format_code) {

      found_string = true;
      break;

    }

    ++offset;

  }

  if (found_string) {

    return offset;

  } else {

    return std::nullopt;

  }

}


