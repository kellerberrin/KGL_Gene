//
// Created by kellerberrin on 27/02/18.
//

#include "kgl_variant_factory.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_genome_vcf_impl.h"
#include "kgl_filter.h"


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
    addGenome(genome_single_variants_, pop_variant_ptr_, variant_quality_);

  }


  reader_ptr_->readVCFFile();


  ExecEnv::log().info("VCF file Records; Read: {}, Rejected: {} (quality={}), Ignored: {} (no matching contig), Error: {}",
                      vcf_record_count_, vcf_record_rejected_, variant_quality_, vcf_record_ignored_, vcf_record_error_);

  ExecEnv::log().info("VCF file Variants; Total generated: {}, Variant database contains :{}, Identical variants ignored: {}",
                      vcf_variant_count_, genome_single_variants_->size(), vcf_variant_count_ - genome_single_variants_->size());

  addGenome(genome_single_variants_, pop_variant_ptr_, variant_quality_);

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
                           genome_single_variants_,
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


void kgl::ParseGenomeVCFImpl::addGenome(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                        std::shared_ptr<PopulationVariant> pop_variant_ptr,
                                        Phred_t read_quality) const {

// Aggregate to compound variants.
  std::shared_ptr<const kgl::GenomeVariant> aggregate_ptr = kgl::VariantFactory::aggregateVariants(genome_db_ptr_,
                                                                                                   genome_name_,
                                                                                                   genome_variant_ptr);

// Filter on  quality >= 5.
  read_quality = read_quality < 5 ? 5 : read_quality;

  std::shared_ptr<const kgl::GenomeVariant> filter_ptr = aggregate_ptr->filterVariants(kgl::QualityFilter(read_quality));

  kgl::ExecEnv::log().info("Filtered for quality: {}, Genome: {} has: {} variants",
                           read_quality, genome_variant_ptr->genomeId(), filter_ptr->size());

// Store the organism variants in the population object.
  pop_variant_ptr->addGenomeVariant(filter_ptr);

}
