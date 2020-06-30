//
// Created by kellerberrin on 23/04/18.
//


#include "kgl_variant_db_unphased.h"
#include "kgl_variant.h"
#include "kel_patterns.h"


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds unphased variants.
// This object hold variants for a contig/chromosome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Use this to copy the object.
std::shared_ptr<kgl::UnphasedContig> kgl::UnphasedContig::deepCopy() const {

  std::shared_ptr<UnphasedContig> contig_copy(std::make_shared<UnphasedContig>(contigId()));

  for(auto const& [offset, variant_vector] : getMap()) {

    for (auto const& variant_count : variant_vector) {

      if (not contig_copy->addVariant(variant_count)) {

        ExecEnv::log().error("UnphasedContig::deepCopy; Cannot add Variant to Contig Copy : {}, at Offset: {}", contig_copy->contigId(), offset);

      }

    }

  }

  return contig_copy;

}

// Unconditionally adds a variant to the contig. Note variants are unordered at a particular offset.
bool kgl::UnphasedContig::addVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  auto result = contig_offset_map_.find(variant_ptr->offset());

  if (result != contig_offset_map_.end()) {
  // Variant offset exists.

    result->second.push_back(variant_ptr);

  } else {
    // add the new offset.
    UnphasedVectorVariantCount offset_vector{variant_ptr};
    auto insert_result = contig_offset_map_.try_emplace(variant_ptr->offset(), offset_vector);

    if (not insert_result.second) {

      ExecEnv::log().error("UnphasedContig::addVariant(); Could not add variant offset: {} to the genome", variant_ptr->offset());
      return false;

    }

  }

  return true;

}


// Test if an equivalent variant already exists in the contig.
bool kgl::UnphasedContig::variantExists(const std::shared_ptr<const Variant>& variant) {

  auto result = contig_offset_map_.find(variant->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.

    for (auto const& existingvariant : result->second) {

      if (existingvariant->equivalent(*variant)) {

        // found the variant.
        return true;

      }
    }
    // If we fall through the loop then no equivalent variant was found.
    return false;

  } else {

    // No variants found at the offset.
    return false;

  }

}

// Adds a unique variant to the contig. No addition if not unique.
bool kgl::UnphasedContig::addUniqueVariant(const std::shared_ptr<const Variant>& variant) {

  if (variantExists(variant)) {
// don't add.

    // Operation normal but variant not unique.
    return true;

  }

  return addVariant(variant);

}


// Counts the variants in a contug.
size_t kgl::UnphasedContig::variantCount() const {


  size_t variant_count = 0;

  for (auto const& offset_variant_vector : getMap()) {

      variant_count += offset_variant_vector.second.size();

  }

  return variant_count;

}


// Creates a copy of the contig that only contains variants passing the filter condition.
std::shared_ptr<kgl::UnphasedContig> kgl::UnphasedContig::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::UnphasedContig> filtered_contig_ptr(std::make_shared<kgl::UnphasedContig>(contigId()));

  for (auto const& [offset, variant_vector] : getMap()) {

    for(auto const& variant_ptr : variant_vector) {

      if (filter.applyFilter(*variant_ptr)) {

        if (not filtered_contig_ptr->addVariant(variant_ptr)) {

          ExecEnv::log().error("UnphasedContig::filterVariants; Problem adding variant at offset: {}, to contig: {}", offset, contigId());

        }

      }

    }

  }

  return filtered_contig_ptr;

}


std::pair<size_t, size_t> kgl::UnphasedContig::validate(const std::shared_ptr<const ContigReference>& contig_db_ptr) const {

  std::pair<size_t, size_t> contig_count{0, 0};

  std::shared_ptr<const DNA5SequenceContig> contig_sequence_ptr = contig_db_ptr->sequence_ptr();

  for (auto const& [offset, variant_vector] : getMap()) {

    contig_count.first += variant_vector.size();

    if (offset >= contig_sequence_ptr->length()) {

      ExecEnv::log().error("UnphasedContig::validate,  Variant offset: {} exceeds total contig: {} size: {}", offset, contig_db_ptr->contigId(), contig_sequence_ptr->length());
      continue;

    }

    for (auto const& variant_ptr : variant_vector) {

      if (not variant_ptr) {

        ExecEnv::log().error("UnphasedContig::validate, Unknown variant: {}", variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));
        continue;

      }

      if (contig_sequence_ptr->subSequence(variant_ptr->offset(), variant_ptr->reference().length()) == variant_ptr->reference()) {

        ++contig_count.second;

      } else {

        ExecEnv::log().error("UnphasedContig::validate, Mismatch, at Contig Offset: {} Sequence is: {}, Variant Reference Sequence is: {}",
                             variant_ptr->offset(),
                             contig_sequence_ptr->subSequence(variant_ptr->offset(), variant_ptr->reference().length()).getSequenceAsString(),
                             variant_ptr->reference().getSequenceAsString());

      }

    }

  }

  return contig_count;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

// Use this to copy the object.
std::shared_ptr<kgl::UnphasedGenome> kgl::UnphasedGenome::deepCopy() const {

  std::shared_ptr<UnphasedGenome> genome_copy(std::make_shared<UnphasedGenome>(genomeId()));

  for (auto const& [contig_id, contig_ptr] :  getMap()) {

    if (not genome_copy->addContig(contig_ptr->deepCopy())) {

      ExecEnv::log().critical("UnphasedGenome::deepCopy(), Genome: {}, Unable to deepcopy Contig: {}", genomeId(), contig_id);

    }

  }

  return genome_copy;

}


// unconditionally merge (retains duplicates) genomes and variants into this genome.
size_t kgl::UnphasedGenome::mergeGenome(const std::shared_ptr<const UnphasedGenome>& merge_genome) {

  if (not merge_genome->processAll(*this, &UnphasedGenome::addVariant)) {

    ExecEnv::log().error("UnphasedGenome::mergeGenome, problem merging genome: {} into genome: {}",
                         merge_genome->genomeId(), genomeId());

  }

  return variantCount();

}


std::shared_ptr<kgl::UnphasedGenome> kgl::UnphasedGenome::uniqueGenome() const {

  std::shared_ptr<UnphasedGenome> unique_genome(std::make_shared<UnphasedGenome>(genomeId()));

  if (not processAll(*unique_genome, &UnphasedGenome::addUniqueVariant)) {

    ExecEnv::log().error("UnphasedGenome::mergeGenome, problem creating unique variant genome from genome: {}", genomeId());

  }

  return unique_genome;

}


bool kgl::UnphasedGenome::addVariant(const std::shared_ptr<const Variant>& variant) {

  std::optional<std::shared_ptr<UnphasedContig>> contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("UnphasedGenome::addVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  if (not contig_opt.value()->addVariant(variant)) {

    ExecEnv::log().error("UnphasedGenome::addVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}


// The first bool is normal operation. The second bool is if a unique variant was added to the genome.
bool kgl::UnphasedGenome::addUniqueVariant(const std::shared_ptr<const Variant>& variant) {

  std::optional<std::shared_ptr<UnphasedContig>> contig_opt = getCreateContig(variant->contigId());
  if (not contig_opt) {

    ExecEnv::log().error("UnphasedGenome::addVariant(), Genome: {} could not get or create Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  bool result = contig_opt.value()->addUniqueVariant(variant);
  if (not result) {

    ExecEnv::log().error("UnphasedGenome::addVariant(), Genome: {} could not add variant to Contig: {}", genomeId(), variant->contigId());
    return false;

  }

  return true;

}



std::optional<std::shared_ptr<kgl::UnphasedContig>> kgl::UnphasedGenome::getCreateContig(const ContigId_t& contig_id) {

  auto result = contig_map_.find(contig_id);

  if (result != contig_map_.end()) {

    return result->second;

  } else {

    std::shared_ptr<UnphasedContig> contig_ptr = std::make_shared<UnphasedContig>(contig_id);
    auto insert_result = contig_map_.try_emplace(contig_id, contig_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("UnphasedGenome::getCreateContig(), Could not add contig: {} to genome : {}", contig_id, genomeId());
      return std::nullopt;

    }

    return contig_ptr;

  }

}


bool kgl::UnphasedGenome::addContig(const std::shared_ptr<UnphasedContig>& contig_ptr) {

  std::pair<ContigId_t, std::shared_ptr<UnphasedContig>> add_contig(contig_ptr->contigId(), contig_ptr);
  auto result = contig_map_.emplace(add_contig);

  if (not result.second) {

    ExecEnv::log().error("UnphasedGenome::addContig(); could not add contig: {} to the genome", contig_ptr->contigId());

  }

  return result.second;

}



size_t kgl::UnphasedGenome::variantCount() const {


  size_t variant_count = 0;

  for (auto contig : getMap()) {

    variant_count += contig.second->variantCount();

  }

  return variant_count;

}



std::shared_ptr<kgl::UnphasedGenome> kgl::UnphasedGenome::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::UnphasedGenome> filtered_genome_ptr(std::make_shared<kgl::UnphasedGenome>(genomeId()));

  for (const auto& contig_variant : getMap()) {

    std::shared_ptr<kgl::UnphasedContig> filtered_contig = contig_variant.second->filterVariants(filter);
    if (not filtered_genome_ptr->addContig(filtered_contig)) {

      ExecEnv::log().critical("UnphasedGenome::filterVariants(), Genome: {}, Unable to inserted filtered Contig: {}", genomeId(), filtered_contig->contigId());

    }

  }

  return filtered_genome_ptr;

}


// Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
// The second integer is the number variants that pass inspection by comparison to the genome database.
std::pair<size_t, size_t> kgl::UnphasedGenome::validate(const std::shared_ptr<const GenomeReference>& genome_db_ptr) const {

  std::pair<size_t, size_t> genome_count{0, 0};
  for (auto const& [contig_id, contig_ptr] : getMap()) {

    std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig_id);

    if (not contig_opt) {

      ExecEnv::log().error("UnphasedGenome::validate, No matching contig found in GenomeDatabase for Variant Contig: {}", contig_id);
      continue;

    }

    std::pair<size_t, size_t> contig_count = contig_ptr->validate(contig_opt.value());

    if (contig_count.first != contig_count.second) {

      ExecEnv::log().warn("UnphasedGenome::validate(), Genome: {} Validation Failed in Contig: {}, Total Variants: {} Validated: {}",
                     genomeId(), contig_id, contig_count.first, contig_count.second);

    }

    genome_count.first += contig_count.first;
    genome_count.second += contig_count.second;

  }

  return genome_count;

}



