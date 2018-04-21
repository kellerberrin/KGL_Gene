//
// Created by kellerberrin on 2/04/18.
//

#include "kgl_variant_factory_vcf_phasing.h"


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
// This object hold variants for a contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VCFContig::addVariant(std::shared_ptr<Variant> variant) {

  auto result = contig_offset_map_.find(variant->offset());

  if (result != contig_offset_map_.end()) {

    result->second.push_back(variant);

    return true;

  } else {

    std::pair<ContigOffset_t, std::vector<std::shared_ptr<Variant>>> new_offset;
    new_offset.first = variant->offset();
    new_offset.second.push_back(variant);
    auto result = contig_offset_map_.insert(new_offset);

    if (not result.second) {

      ExecEnv::log().error("VCFContig::addVariant(), Could not add variant offset: {} to the genome", variant->offset());
      return false;

    }

    return true;

  }

}


size_t kgl::VCFContig::variantCount() const {


  size_t variant_count = 0;

  for (auto offset : getMap()) {

    variant_count += offset.second.size();

  }

  return variant_count;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
// This object hold variants for a genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VCFGenome::addVariant(std::shared_ptr<Variant> variant) {

  std::shared_ptr<VCFContig> contig_ptr;
  getCreateContig(variant->contigId(), contig_ptr);

  contig_ptr->addVariant(variant);

  return true;

}


bool kgl::VCFGenome::getCreateContig(const ContigId_t& contig_id, std::shared_ptr<VCFContig>& contig_ptr) {

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    contig_ptr = result->second;
    return true;

  } else {

    contig_ptr = std::make_shared<VCFContig>(contig_id);
    std::pair<ContigId_t, std::shared_ptr<VCFContig>> new_contig(contig_id, contig_ptr);
    auto result = contig_map_.insert(new_contig);

    if (not result.second) {

      ExecEnv::log().critical("VCFGenome::getCreateContig(), Serious Error, could not add contig: {} to the genome", contig_id);

    }

    return result.second;

  }

}




size_t kgl::VCFGenome::variantCount() const {


  size_t variant_count = 0;

  for (auto contig : getMap()) {

    variant_count += contig.second->variantCount();

  }

  return variant_count;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds multi ploid variants until they can be phased.
// This object hold variants for a population.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VCFPopulation::getCreateGenome(const GenomeId_t& genome_id,
                                         std::shared_ptr<VCFGenome>& genome) {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    genome = result->second;
    return true;

  } else {

    genome = std::make_shared<VCFGenome>(genome_id);
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
bool kgl::GenomePhasing::haploidPhasing(std::shared_ptr<const VCFPopulation> vcf_population_ptr,
                                        std::shared_ptr<const GenomeDatabase> genome_db,
                                        std::shared_ptr<PopulationVariant> haploid_population)  {

  bool result = true;

  for (auto genome : vcf_population_ptr->getMap()) {

    // Create the GenomeVariant object.
    std::shared_ptr<GenomeVariant> genome_variant = GenomeVariant::emptyGenomeVariant(genome.first, genome_db);

    // Add all the contig.
    for (auto contig : genome.second->getMap()) {

      // Add all offsets.
      for (auto offset : contig.second->getMap()) {

        // add all the variants
        for (auto variant : offset.second) {

          genome_variant->addVariant(variant);

        } // All variants.

      } // All offsets

    } // All contigs

    if (not haploid_population->addGenomeVariant(genome_variant)) {

      ExecEnv::log().error("Unable to add genome: {} to haploid population", genome_variant->genomeId());
      result = false;

    }

  } // All genomes

  ExecEnv::log().info("Haploid Phasing, pre_phase variants: {}, genomes: {}, resultant population variants: {}, genomes: {}",
                      vcf_population_ptr->variantCount(), vcf_population_ptr->getMap().size(),
                      haploid_population->variantCount(), haploid_population->getMap().size());

  return result;

}
