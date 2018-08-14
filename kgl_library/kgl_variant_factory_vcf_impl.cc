//
// Created by kellerberrin on 22/01/18.
//


#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


void kgl::ParseVCFImpl::readParseVCFImpl() {


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

  }

  // multi-threaded
  reader_ptr_->readVCFFile();
  // single threaded

  ExecEnv::log().info("VCF file Records; Read: {}, Rejected: {} (quality={}), Ignored: {} (no matching contig), Error: {}",
                      vcf_record_count_, vcf_record_rejected_, variant_quality_, vcf_record_ignored_, vcf_record_error_);



}



size_t kgl::ParseVCFImpl::addThreadSafeGenomeVariant(std::shared_ptr<const Variant>& variant_ptr) {

  AutoMutex auto_mutex(mutex_);

  std::shared_ptr<UnphasedGenome> genome;
  unphased_population_ptr_->getCreateGenome(variant_ptr->genomeId(), genome); // thread safe
  genome->addVariant(variant_ptr); // thread safe

  return 1;

}


// Set up the genomes first rather than on-the-fly.
// Some genomes may have no variants (the model/reference genome) and thus would not be created.
void kgl::ParseVCFImpl::setupVCFPopulation() {

  AutoMutex auto_mutex(mutex_);

  for (auto genome_id : getGenomeNames())  {

    std::shared_ptr<UnphasedGenome> genome;
    unphased_population_ptr_->getCreateGenome(genome_id, genome);

  }

}

