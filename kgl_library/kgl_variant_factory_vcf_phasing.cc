//
// Created by kellerberrin on 2/04/18.
//

#include "kgl_variant_factory_vcf_phasing.h"


namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Thread safe. This object can be used by multiple consumer threads.
// An internal parser variant object that holds multi ploid variants until they can be phased.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VCFGenome::addVariant(std::shared_ptr<Variant> variant) {

  std::shared_ptr<ContigMultiMap> contig_ptr;
  getCreateContig(variant->contigId(), contig_ptr);

  contig_ptr->insert(ContigMultiMap::value_type(variant->offset(), variant));

  return true;

}


bool kgl::VCFGenome::getCreateContig(const ContigId_t& contig_id, std::shared_ptr<ContigMultiMap>& contig_ptr) {

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    contig_ptr = result->second;
    return true;

  } else {

    contig_ptr = std::make_shared<ContigMultiMap>();
    std::pair<ContigId_t, std::shared_ptr<ContigMultiMap>> new_contig(contig_id, contig_ptr);
    auto result = contig_map_.insert(new_contig);

    if (not result.second) {

      ExecEnv::log().critical("VCFGenome::getCreateContig(), Serious Error, could not add contig: {} to the genome", contig_id);

    }

    return result.second;

  }

}




size_t kgl::VCFGenome::variantCount() const {


  size_t variant_count = 0;

  for (auto contig : contig_map_) {

    variant_count += contig.second->size();

  }

  return variant_count;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Thread safe. This object can be used by multiple consumer threads.
// An internal parser variant object that holds multi ploid variants until they can be phased.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VCFPopulation::getCreateGenome(const GenomeId_t& genome_id,
                                         std::shared_ptr<VCFGenome>& genome) {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    genome = result->second;
    return true;

  } else {

    genome = std::make_shared<VCFGenome>();
    std::pair<GenomeId_t, std::shared_ptr<VCFGenome>> new_genome(genome_id, genome);
    auto result = genome_map_.insert(new_genome);

    if (not result.second) {

      ExecEnv::log().critical("VCFPopulation::getCreateGenome(), Serious Error, could not add genome: {} to the population", genome_id);

    }

    return result.second;

  }

}


size_t kgl::VCFPopulation::variantCount() const {

  size_t variant_count = 0;

  for (auto genome : genome_map_) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Not Thread safe.
// This object accepts holds multi-ploid variants and phases them.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


// Just copy into a population object.
bool kgl::GenomePhasing::haploidPhasing(const VCFPopulation& vcf_population,
                                        std::shared_ptr<const GenomeDatabase> genome_db,
                                        std::shared_ptr<PopulationVariant> haploid_population)  {

  bool result = true;

  for (auto genome : vcf_population.getMap()) {

    // Create the GenomeVariant object.
    std::shared_ptr<GenomeVariant> genome_variant = GenomeVariant::emptyGenomeVariant(genome.first, genome_db);

    // Add all the contig.
    for (auto contig : genome.second->getMap()) {

      // Add all variants.
      for (auto variant : *(contig.second)) {

        genome_variant->addVariant(variant.second);

      }

    }

    if (not haploid_population->addGenomeVariant(genome_variant)) {

      ExecEnv::log().error("Unable to add genome: {} to haploid population", genome_variant->genomeId());
      result = false;

    }

  }

  ExecEnv::log().info("Haploid Phasing, pre_phase variants: {}, genomes: {}, resultant population variants: {}, genomes: {}",
                      vcf_population.variantCount(), vcf_population.getMap().size(),
                      haploid_population->variantCount(), haploid_population->getMap().size());

  return result;

}
