//
// Created by kellerberrin on 22/01/18.
//


#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_vcf.h"


namespace kgl = kellerberrin::genome;



size_t kgl::ParseVCFImpl::addThreadSafeGenomeVariant(std::shared_ptr<const Variant>& variant_ptr) {

  AutoMutex auto_mutex(mutex_); // Write Locked

  std::shared_ptr<UnphasedGenome> genome;

  if (not unphased_population_ptr_->getCreateGenome(variant_ptr->genomeId(), genome)) {

    ExecEnv::log().error("ParseVCFImpl::addThreadSafeGenomeVariant; Could not add/create genome: {}", variant_ptr->genomeId());
    return 0;

  }

  if (not genome->addVariant(variant_ptr)) { // thread safe

    ExecEnv::log().error("ParseVCFImpl::addThreadSafeGenomeVariant; Could not add variant to genome: {}", variant_ptr->genomeId());
    return 0;

  }

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

    if (not unphased_population_ptr_->getCreateGenome(genome_id, genome_ptr)) {

      ExecEnv::log().critical("Could not create genome: {} in the unphased population", genome_id);

    }

    for (auto contig : genome_db_ptr->getMap()) {

      std::shared_ptr<UnphasedContig> contig_ptr = nullptr;

      if (not genome_ptr->getCreateContig(contig.first, contig_ptr)) {

        ExecEnv::log().critical("Could not create contig: {} in genome: {} in the unphased population", contig.first, genome_id);

      }

    }

  }

}

bool kgl::ParseVCFImpl::parseCigarItems(const std::string& genome_name,
                                      std::shared_ptr<const ContigFeatures> contig_ptr,
                                      const std::vector<CigarEditItem>&, // parsed_cigar,  // Cigar not used
                                      ContigOffset_t contig_offset,
                                      const std::string& reference_text,
                                      const std::string& alternate_text,
                                      std::shared_ptr<const VariantEvidence> evidence_ptr,
                                      size_t& record_variants)  {

  StringDNA5 reference_str(reference_text);
  StringDNA5 alternate_str(alternate_text);

  std::shared_ptr<const Variant> variant_ptr(std::make_shared<VCFVariant>(genome_name,
                                                                          contig_ptr->contigId(),
                                                                          VariantSequence::UNPHASED,
                                                                          contig_offset,
                                                                          evidence_ptr,
                                                                          std::move(reference_str),
                                                                          std::move(alternate_str)));

  record_variants = addThreadSafeGenomeVariant(variant_ptr);

  return true;

}

