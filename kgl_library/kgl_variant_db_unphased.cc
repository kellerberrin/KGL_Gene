//
// Created by kellerberrin on 23/04/18.
//


#include "kgl_variant_db_unphased.h"


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::UnphasedContig::addVariant(std::shared_ptr<Variant> variant) {

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

      ExecEnv::log().error("UnphasedContig::addVariant(), Could not add variant offset: {} to the genome", variant->offset());
      return false;

    }

    return true;

  }

}


size_t kgl::UnphasedContig::variantCount() const {


  size_t variant_count = 0;

  for (auto offset : getMap()) {

    variant_count += offset.second.size();

  }

  return variant_count;

}

// true if vector.size() < 2
bool kgl::UnphasedContig::isHomozygous(const std::vector<std::shared_ptr<Variant>>& variant_vector) {

  bool is_homozygous = true;
  for (auto variant : variant_vector) {

    if (not variant_vector.front()->equivalent(*variant)) {

      is_homozygous = false;
      break;

    }

  }

  return is_homozygous;

}


bool kgl::UnphasedContig::removeConflictingVariants() {

  for (auto& variant_vector : contig_offset_map_) {

    if (not variant_vector.second.empty()) {
      // remove all but the first SNP variant
      auto variant_iter = variant_vector.second.begin();
      while (variant_iter != variant_vector.second.end()) {

        if ((*variant_iter)->isSNP()) break;
        ++variant_iter;

      }

      if (variant_iter != variant_vector.second.end()) {

        std::shared_ptr<Variant> first_variant = *variant_iter;
        variant_vector.second.clear();
        variant_vector.second.push_back(first_variant); // add the first variant back into the vector.

      } else {

        contig_offset_map_.erase(variant_vector.first); // erase the vector.

      }


    } else {

      ExecEnv::log().warn("removeConflictingVariants(); empty variant vector in unphased variant contig: {}, offset: {}",
                          contig_id_, variant_vector.first);
      contig_offset_map_.erase(variant_vector.first); // erase the empty vector.

    }

  }

  return true;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::UnphasedGenome::addVariant(std::shared_ptr<Variant> variant) {

  std::shared_ptr<UnphasedContig> contig_ptr;
  getCreateContig(variant->contigId(), contig_ptr);

  contig_ptr->addVariant(variant);

  return true;

}


bool kgl::UnphasedGenome::getCreateContig(const ContigId_t& contig_id, std::shared_ptr<UnphasedContig>& contig_ptr) {

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    contig_ptr = result->second;
    return true;

  } else {

    contig_ptr = std::make_shared<UnphasedContig>(contig_id);
    std::pair<ContigId_t, std::shared_ptr<UnphasedContig>> new_contig(contig_id, contig_ptr);
    auto result = contig_map_.insert(new_contig);

    if (not result.second) {

      ExecEnv::log().critical("UnphasedGenome::getCreateContig(), Serious Error, could not add contig: {} to the genome", contig_id);

    }

    return result.second;

  }

}

bool kgl::UnphasedGenome::removeConflictingVariants() {

  for (auto contig : getMap()) {

    contig.second->removeConflictingVariants();

  }

  return true;

}


size_t kgl::UnphasedGenome::variantCount() const {


  size_t variant_count = 0;

  for (auto contig : getMap()) {

    variant_count += contig.second->variantCount();

  }

  return variant_count;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a population.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::UnphasedPopulation::getCreateGenome(const GenomeId_t& genome_id,
                                              std::shared_ptr<UnphasedGenome>& genome) {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    genome = result->second;
    return true;

  } else {

    genome = std::make_shared<UnphasedGenome>(genome_id);
    std::pair<GenomeId_t, std::shared_ptr<UnphasedGenome>> new_genome(genome_id, genome);
    auto result = genome_map_.insert(new_genome);

    if (not result.second) {

      ExecEnv::log().critical("UnphasedPopulation::getCreateGenome(), Serious Error, could not add genome: {} to the population", genome_id);

    }

    return result.second;

  }

}


bool kgl::UnphasedPopulation::removeConflictingVariants() {

  for (auto genome : getMap()) {

    genome.second->removeConflictingVariants();

  }

  return true;

}


size_t kgl::UnphasedPopulation::variantCount() const {

  size_t variant_count = 0;

  for (auto genome : genome_map_) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}


bool kgl::UnphasedPopulation::genomePhasingStats(const GenomeId_t& genome_id,
                                                 bool snp_only,
                                                 size_t& heterozygous,
                                                 size_t& homozygous,
                                                 size_t& singleheterozygous) const {

  bool return_result = true;

  heterozygous = 0;
  homozygous = 0;
  singleheterozygous = 0;

  auto result = genome_map_.find(genome_id);

  if (result == genome_map_.end()) {

    ExecEnv::log().error("UnphasedPopulation::genomePhasingStats(); Could not find genome: {}", genome_id);
    return false;

  }

  for (auto contig : result->second->getMap()) {

    for (auto offset : contig.second->getMap()) {

      if (offset.second.size() == 0) {

        ExecEnv::log().error("UnphasedPopulation::genomePhasingStats(); Zero sized variant vector, genome: {}, contig: {}, offset: {}",
                             genome_id, contig.first, offset.first);
        return_result = false;
        continue;

      }

      if (snp_only and not offset.second.front()->isSNP()) {

        continue;

      }

      bool is_homozygous = UnphasedContig::isHomozygous(offset.second);

      if (offset.second.size() >= 2 and not is_homozygous) {

        ++heterozygous;

      }
      if (offset.second.size() >= 2 and is_homozygous) {

        ++homozygous;

      }
      if (offset.second.size() == 1) {

        ++singleheterozygous;

      }

    } // offset

  } // contig

  return return_result;

}


bool kgl::UnphasedPopulation::getUnphasedVariants(const GenomeId_t& genome_id,
                                                  const ContigId_t& contig_id,
                                                  ContigOffset_t offset,
                                                  std::vector<std::shared_ptr<Variant>>& variant_vector) const {

  auto genome_result = genome_map_.find(genome_id);

  if (genome_result == genome_map_.end()) {

    ExecEnv::log().error("UnphasedPopulation::getUnphasedVariants(); Could not find genome: {}", genome_id);
    return false;

  }

  auto contig_result = genome_result->second->getMap().find(contig_id);

  if (contig_result == genome_result->second->getMap().end()) {

    ExecEnv::log().error("UnphasedPopulation::getUnphasedVariants(); Could not find contig: {}", contig_id);
    return false;

  }

  auto offset_result = contig_result->second->getMap().find(offset);

  if (offset_result == contig_result->second->getMap().end()) {

    ExecEnv::log().error("UnphasedPopulation::getUnphasedVariants(); Could not find offset: {}", offset);
    return false;

  }

  variant_vector = offset_result->second;

  return true;

}


void kgl::UnphasedPopulation::popStatistics() const {

  size_t total_variants = 0;

  for (auto genome : getMap()) {

    size_t genome_variant_count = genome.second->variantCount();

    ExecEnv::log().info("Genome: {}, Unphased Variant Count:{}", genome.first, genome_variant_count);

    total_variants += genome_variant_count;

  }

  ExecEnv::log().info("Total Unphased Variant Count:{}", total_variants);

}

std::vector<kgl::GenomeId_t> kgl::UnphasedPopulation::genomeList() const {

  std::vector<kgl::GenomeId_t> genome_list;

  for (auto genome : getMap()) {

    genome_list.push_back(genome.first);

  }

  return genome_list;

}
