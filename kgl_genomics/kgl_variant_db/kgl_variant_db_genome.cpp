//
// Created by kellerberrin on 25/12/20.
//

#include "kgl_variant_db_genome.h"
#include "kgl_variant_filter_db_variant.h"

namespace kgl = kellerberrin::genome;


// Use this to copy the object.
std::shared_ptr<kgl::GenomeDB> kgl::GenomeDB::deepCopy() const {

  // Can use the shallow filter because all variants are copied across.
  return viewFilter(TrueFilter());

}

// unconditionally merge (retains duplicates) genomes and variants into this genome.
size_t kgl::GenomeDB::mergeGenome(const std::shared_ptr<const GenomeDB>& merge_genome) {

  if (not merge_genome->processAll(*this, &GenomeDB::addVariant)) {

    ExecEnv::log().error("GenomeDB::mergeGenome, problem merging genome: {} into genome: {}",
                         merge_genome->genomeId(), genomeId());

  }

  return variantCount();

}


bool kgl::GenomeDB::addVariant(const std::shared_ptr<const Variant>& variant) {

  auto contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("GenomeDB::addVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  if (not contig_opt.value()->addVariant(variant)) {

    ExecEnv::log().error("GenomeDB::addVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}


std::optional<std::shared_ptr<kgl::ContigDB>> kgl::GenomeDB::getCreateContig(const ContigId_t& contig_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto found = contig_map_.find(contig_id);

  if (found != contig_map_.end()) {

    return found->second;

  } else {

    auto contig_ptr = std::make_shared<ContigDB>(contig_id);
    auto [inserted, insert_result] = contig_map_.insert({contig_id, contig_ptr});

    if (not insert_result) {

      ExecEnv::log().error("GenomeDB::getCreateContig(), Could not add contig_ref_ptr: {} to genome : {}", contig_id, genomeId());
      return std::nullopt;

    }

    return inserted->second;

  }

}


std::optional<std::shared_ptr<kgl::ContigDB>> kgl::GenomeDB::getContig(const ContigId_t& contig_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}

std::optional<std::shared_ptr<const kgl::ContigDB>> kgl::GenomeDB::getContig(const ContigId_t& contig_id) const {

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


bool kgl::GenomeDB::addContig(std::shared_ptr<ContigDB> contig_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto [iter, result] = contig_map_.try_emplace(contig_ptr->contigId(), contig_ptr);

  if (not result) {

    ExecEnv::log().error("GenomeDB::addContig(); could not add contig_ref_ptr: {} (duplicate) to the genome", contig_ptr->contigId());

  }

  return result;

}


size_t kgl::GenomeDB::variantCount() const {

  size_t variant_count{0};

  for (auto const& [contig_id, contig_ptr] : getMap()) {

    variant_count += contig_ptr->variantCount();

  }

  return variant_count;

}


std::unique_ptr<kgl::GenomeDB> kgl::GenomeDB::viewFilter(const BaseFilter& filter) const {

  // Only genome filter is implemented at this level.
  if (filter.filterType() == FilterBaseType::GENOME_FILTER) {

    const FilterGenomes& genome_filter = static_cast<const FilterGenomes&>(filter);
    return genome_filter.applyFilter(*this);

  }
  // All other filters.
  // Filter the contigs.
  std::unique_ptr<GenomeDB> filtered_genome_ptr(std::make_unique<GenomeDB>(genomeId()));
  for (const auto& [contig_id, contig_ptr] : getMap()) {

    std::shared_ptr<ContigDB> filtered_contig = contig_ptr->viewFilter(filter);
    if (not filtered_genome_ptr->addContig(filtered_contig)) {

      ExecEnv::log().error("GenomeDB::filter; Genome: {}, Unable to insert filtered Contig: {}", genomeId(), contig_id);

    }

  }

  return filtered_genome_ptr;

}


// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
std::pair<size_t, size_t> kgl::GenomeDB::selfFilter(const BaseFilter& filter) {

  // Only genome filter is implemented at this level.
  if (filter.filterType() == FilterBaseType::GENOME_FILTER) {

    size_t prior_count = variantCount();

    std::shared_ptr<const FilterGenomes> genome_filter = std::dynamic_pointer_cast<const FilterGenomes>(filter.clone());
    auto genome_ptr = genome_filter->applyFilter(*this);
    contig_map_ = std::move(genome_ptr->contig_map_);

    size_t post_count = variantCount();

    return {prior_count, post_count};

  }

  // All other filters.
  // Filter the contigs.
  std::pair<size_t, size_t> genome_count{0, 0};
  for (auto& [contig_id, contig_ptr] : contig_map_) {

    auto contig_count = contig_ptr->selfFilter(filter);
    genome_count.first += contig_count.first;
    genome_count.second += contig_count.second;

  }

  return genome_count;

}

// Deletes any empty Contigs, returns number deleted.
size_t kgl::GenomeDB::trimEmpty() {

  size_t delete_count{0};
  // Delete empty genomes.
  auto it = contig_map_.begin();
  while (it != contig_map_.end()) {

    if (it->second->variantCount() == 0) {

      it = contig_map_.erase(it);
      ++delete_count;

    } else {

      ++it;

    }

  }

  return delete_count;

}

// Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
// The second integer is the number variants that pass inspection by comparison to the genome database.
std::pair<size_t, size_t> kgl::GenomeDB::validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const {

  std::pair<size_t, size_t> genome_count{0, 0};
  for (auto const& [contig_id, contig_ptr] : getMap()) {

    std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig_id);

    if (not contig_opt) {

      ExecEnv::log().error("GenomeDB::validate, No matching contig_ref_ptr found in GenomeDatabase for Variant Contig: {}", contig_id);
      continue;

    }

    std::pair<size_t, size_t> contig_count = contig_ptr->validate(contig_opt.value());

    if (contig_count.first != contig_count.second) {

      ExecEnv::log().warn("GenomeDB::validate(), Genome: {} Validation Failed in Contig: {}, Total filter: {} Validated: {}",
                          genomeId(), contig_id, contig_count.first, contig_count.second);

    }

    genome_count.first += contig_count.first;
    genome_count.second += contig_count.second;

  }

  return genome_count;

}


// Create an equivalent genome that has canonical variants, SNP are represented by '1X', Deletes by '1MnD'
// and Inserts by '1MnI'. The genome structure is re-created and is not a shallow copy.
std::unique_ptr<kgl::GenomeDB> kgl::GenomeDB::canonicalGenome() const {

  // Create the new genome.
  std::unique_ptr<GenomeDB> canonical_genome_ptr(std::make_unique<GenomeDB>(genomeId()));
  // Populate with canonical contigs.
  for (auto const& [contig_id, contig_ptr] : getMap()) {

    std::shared_ptr<ContigDB> canonical_contig_ptr = contig_ptr->canonicalContig();
    if (not canonical_genome_ptr->addContig(canonical_contig_ptr)) {

      ExecEnv::log().error("GenomeDB::canonicalGenome; cannot add contig_ref_ptr: {} to genome: {}", contig_id, genomeId());

    }

  }

  return canonical_genome_ptr;

}


bool kgl::GenomeDB::processAll(const VariantProcessFunc& procFunc) const {


  for (auto const& [contig, contig_ptr] : getMap()) {

    if (not contig_ptr->processAll(procFunc)) {

      ExecEnv::log().error("GenomeDB::processAll; Problem executing general purpose function for contig_ref_ptr: {}", contig);
      return false;

    }

  }

  return true;

}


