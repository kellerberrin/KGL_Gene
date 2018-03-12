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

  // Now in single threaded code so transfer the variants to the population object.
  for (auto genome : thread_safe_population_.getMap()) {

    pop_variant_ptr_->addGenomeVariant(genome.second);

  }

}



size_t kgl::ParseVCFImpl::addThreadSafeGenomeVariant(std::shared_ptr<GenomeVariant> genome_variants,
                                                     std::shared_ptr<const Variant> variant_ptr) const {

  AutoMutex auto_mutex(mutex_);

  genome_variants->addVariant(variant_ptr);

  return 1;

}



bool kgl::ThreadSafePopulation::getCreateGenomeVariant(const GenomeId_t& genome_id,
                                                       const std::shared_ptr<const GenomeDatabase> genome_db,
                                                       std::shared_ptr<GenomeVariant>& genome_variant) {
  AutoMutex auto_mutex(mutex_);

  auto result = population_variant_map_.find(genome_id);

  if (result != population_variant_map_.end()) {

    genome_variant = result->second;
    return true;

  } else {

    genome_variant = GenomeVariant::emptyGenomeVariant(genome_id, genome_db);
    auto result = population_variant_map_.insert(std::pair<GenomeId_t, std::shared_ptr<GenomeVariant>>(genome_id, genome_variant));

    if (not result.second) {

      ExecEnv::log().error("Serious Error, could not add genome: {} to the population", genome_id);
      genome_variant = nullptr;

    }

    return result.second;

  }

}


size_t kgl::ThreadSafePopulation::variantCount() const {

  AutoMutex auto_mutex(mutex_);

  size_t variant_count = 0;

  for (auto genome : population_variant_map_) {

    variant_count += genome.second->size();

  }

  return variant_count;

}
