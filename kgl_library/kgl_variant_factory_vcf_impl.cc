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

  ExecEnv::log().info("VCF file Records; Read: {}, Rejected: {}, Ignored: {} (no matching contig), Error: {}",
                      vcf_record_count_, vcf_record_rejected_, vcf_record_ignored_, vcf_record_error_);



}



size_t kgl::ParseVCFImpl::addThreadSafeGenomeVariant(std::shared_ptr<const Variant>& variant_ptr) {

  AutoMutex auto_mutex(mutex_);

  std::shared_ptr<UnphasedGenome> genome;
  unphased_population_ptr_->getCreateGenome(variant_ptr->genomeId(), genome); // thread safe
  genome->addVariant(variant_ptr); // thread safe

  return 1;

}


// Set up the genomes/contigs first rather than on-the-fly.
// Some genomes may have no variants (e.g. the model/reference genome 3D7)
// and thus these genomes/contigs would not be created on-the-fly.
void kgl::ParseVCFImpl::setupPopulationStructure(std::shared_ptr<const GenomeDatabase> genome_db_ptr) {

  AutoMutex auto_mutex(mutex_);

  ExecEnv::log().info("setupPopulationStructure; Creating a population of {} genomes and {} contigs",
                      getGenomeNames().size(), genome_db_ptr->getMap().size());

  for (auto genome_id : getGenomeNames())  {

    std::shared_ptr<UnphasedGenome> genome_ptr = nullptr;
    unphased_population_ptr_->getCreateGenome(genome_id, genome_ptr);

    if (not genome_ptr) {

      ExecEnv::log().critical("Could not create genome: {} in the unphased population", genome_id);

    }

    for (auto contig : genome_db_ptr->getMap()) {

      std::shared_ptr<UnphasedContig> contig_ptr = nullptr;
      genome_ptr->getCreateContig(contig.first, contig_ptr);

      if (not contig_ptr) {

        ExecEnv::log().critical("Could not create contig: {} in genome: {} in the unphased population", contig.first, genome_id);

      }

    }

  }

}

