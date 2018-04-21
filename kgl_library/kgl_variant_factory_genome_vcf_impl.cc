//
// Created by kellerberrin on 27/02/18.
//

#include "kgl_variant_factory.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_genome_vcf_impl.h"
#include "kgl_filter.h"
#include "kgl_variant_phasing_statistics.h"
#include "kgl_vcf_parser_data.h"


namespace kgl = kellerberrin::genome;



void kgl::ParseGenomeVCFImpl::readParseVCFImpl() {

  // Process records.
  // Copy the file record by record.
  vcf_record_count_ = 0;
  vcf_record_error_ = 0;
  vcf_record_ignored_ = 0;
  vcf_record_rejected_ = 0;
  vcf_variant_count_ = 0;


  // Investigate header.
  ActiveContigMap active_contig_map;
  if (not ParseVCFMiscImpl::parseVcfHeader(genome_db_ptr_, reader_ptr_->readHeader(), active_contig_map, false)) {

    ExecEnv::log().error("Problem parsing header information in VCF file: {}. No variants processed.", vcf_file_name_);
    return;

  }


  reader_ptr_->readVCFFile();


  ExecEnv::log().info("VCF file Records; Read: {}, Rejected: {} (quality={}), Ignored: {} (no matching contig), Error: {}",
                      vcf_record_count_, vcf_record_rejected_, variant_quality_, vcf_record_ignored_, vcf_record_error_);

  ExecEnv::log().info("VCF file Variants; Total generated: {}", vcf_variant_count_);


}


// This is multithreaded code called from the reader defined above.
void kgl::ParseGenomeVCFImpl::ProcessVCFRecord(const seqan::VcfRecord& record)
{

  size_t record_variants = 0;

  ++vcf_record_count_;

  ContigId_t contig_id = reader_ptr_->getContig(record.rID);
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (genome_db_ptr_->getContigSequence(contig_id, contig_ptr)) {

    bool record_quality_ok;

    if (not parseVcfRecord(genome_name_,
                           record,
                           contig_ptr,
                           variant_quality_,
                           record_quality_ok,
                           record_variants)) {


      ++vcf_record_error_;
      ExecEnv::log().error("Error parsing VCF record");

    }

    if (not record_quality_ok) {

      ++vcf_record_rejected_;

    }

  } else {

    ++vcf_record_ignored_;

  }


  for (size_t idx = 0; idx < record_variants; ++idx) {

    ++vcf_variant_count_;

    if (vcf_variant_count_ % VARIANT_REPORT_INTERVAL_ == 0) {

      ExecEnv::log().info("VCF file, generated: {} variants", vcf_variant_count_);

    }

  }

}

