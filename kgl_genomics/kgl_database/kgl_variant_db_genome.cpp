//
// Created by kellerberrin on 25/12/20.
//

#include "kgl_variant_db_genome.h"

namespace kgl = kellerberrin::genome;




// Use this to copy the object.
std::shared_ptr<kgl::GenomeVariantArray> kgl::GenomeVariantArray::deepCopy() const {

  std::shared_ptr<GenomeVariantArray> genome_copy(std::make_shared<GenomeVariantArray>(genomeId()));

  for (auto const& [contig_id, contig_ptr] :  getMap()) {

    if (not genome_copy->addContig(contig_ptr->deepCopy())) {

      ExecEnv::log().critical("GenomeVariantArray::deepCopy(), Genome: {}, Unable to deepcopy Contig: {}", genomeId(), contig_id);

    }

  }

  return genome_copy;

}


// unconditionally merge (retains duplicates) genomes and variants into this genome.
size_t kgl::GenomeVariantArray::mergeGenome(const std::shared_ptr<const GenomeVariantArray>& merge_genome) {

  if (not merge_genome->processAll(*this, &GenomeVariantArray::addVariant)) {

    ExecEnv::log().error("GenomeVariantArray::mergeGenome, problem merging genome: {} into genome: {}",
                         merge_genome->genomeId(), genomeId());

  }

  return variantCount();

}


bool kgl::GenomeVariantArray::addVariant(const std::shared_ptr<const Variant>& variant) {

  auto contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("GenomeVariantArray::addVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  if (not contig_opt.value()->addVariant(variant)) {

    ExecEnv::log().error("GenomeVariantArray::addVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}

bool kgl::GenomeVariantArray::addUniqueUnphasedVariant(const std::shared_ptr<const Variant>& variant) {

  auto contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("GenomeVariantArray::addUniqueUnphasedVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  if (not contig_opt.value()->addUniqueUnphasedVariant(variant)) {

    ExecEnv::log().error("GenomeVariantArray::addUniqueUnphasedVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}


std::optional<std::shared_ptr<kgl::ContigOffsetVariant>> kgl::GenomeVariantArray::getCreateContig(const ContigId_t& contig_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    auto contig_ptr = std::make_shared<ContigOffsetVariant>(contig_id);
    auto insert_result = contig_map_.try_emplace(contig_id, contig_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("GenomeVariantArray::getCreateContig(), Could not add contig: {} to genome : {}", contig_id, genomeId());
      return std::nullopt;

    }

    return contig_ptr;

  }

}


std::optional<std::shared_ptr<kgl::ContigOffsetVariant>> kgl::GenomeVariantArray::getContig(const ContigId_t& contig_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}

std::optional<std::shared_ptr<const kgl::ContigOffsetVariant>> kgl::GenomeVariantArray::getContig(const ContigId_t& contig_id) const {

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


bool kgl::GenomeVariantArray::addContig(const std::shared_ptr<ContigOffsetVariant>& contig_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_map_.try_emplace(contig_ptr->contigId(), contig_ptr);

  if (not result.second) {

    ExecEnv::log().error("GenomeVariantArray::addContig(); could not add contig: {} (duplicate) to the genome", contig_ptr->contigId());

  }

  return result.second;

}


size_t kgl::GenomeVariantArray::variantCount() const {

  size_t variant_count = 0;

  for (auto contig : getMap()) {

    variant_count += contig.second->variantCount();

  }

  return variant_count;

}


std::shared_ptr<kgl::GenomeVariantArray> kgl::GenomeVariantArray::filterVariants(const VariantFilter& filter) const {

  std::shared_ptr<GenomeVariantArray> filtered_genome_ptr(std::make_shared<GenomeVariantArray>(genomeId()));

  for (const auto& [contig_id, contig_ptr] : getMap()) {

    auto filtered_contig = contig_ptr->filterVariants(filter);
    if (not filtered_genome_ptr->addContig(filtered_contig)) {

      ExecEnv::log().error("GenomeVariantArray::filterVariants(), Genome: {}, Unable to inserted filtered Contig: {}", genomeId(), contig_id);

    }

  }

  return filtered_genome_ptr;

}


// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
std::pair<size_t, size_t> kgl::GenomeVariantArray::inSituFilter(const VariantFilter& filter) {

  // Filter the contigs.
  std::pair<size_t, size_t> genome_count{0, 0};
  for (auto& [contig_id, contig_ptr] : contig_map_) {

    auto contig_count = contig_ptr->inSituFilter(filter);
    genome_count.first += contig_count.first;
    genome_count.second += contig_count.second;

  }
  // Delete the empty contigs.
  for (auto it = contig_map_.begin(); it != contig_map_.end(); ++it) {

    if (it->second->getMap().empty()) {

      it = contig_map_.erase(it);

    }

  }

  return genome_count;

}


// Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
// The second integer is the number variants that pass inspection by comparison to the genome database.
std::pair<size_t, size_t> kgl::GenomeVariantArray::validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const {

  std::pair<size_t, size_t> genome_count{0, 0};
  for (auto const& [contig_id, contig_ptr] : getMap()) {

    std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig_id);

    if (not contig_opt) {

      ExecEnv::log().error("GenomeVariantArray::validate, No matching contig found in GenomeDatabase for Variant Contig: {}", contig_id);
      continue;

    }

    std::pair<size_t, size_t> contig_count = contig_ptr->validate(contig_opt.value());

    if (contig_count.first != contig_count.second) {

      ExecEnv::log().warn("GenomeVariantArray::validate(), Genome: {} Validation Failed in Contig: {}, Total Variants: {} Validated: {}",
                          genomeId(), contig_id, contig_count.first, contig_count.second);

    }

    genome_count.first += contig_count.first;
    genome_count.second += contig_count.second;

  }

  return genome_count;

}


bool kgl::GenomeVariantArray::getSortedVariants( ContigId_t contig_id,
                                                 PhaseId_t phase,
                                                 ContigOffset_t start,
                                                 ContigOffset_t end,
                                                 OffsetVariantMap &variant_map) const {


  auto result = contig_map_.find(contig_id);

  if (result == contig_map_.end()) {

    ExecEnv::log().error("Contig Id: {} not found in Genome Variant: {}", contig_id, genomeId());
    return false;

  }

  std::shared_ptr<ContigOffsetVariant> contig_ptr = result->second;

  return contig_ptr->getSortedVariants(phase, start, end, variant_map);

}


