//
// Created by kellerberrin on 25/12/20.
//


#include "kgl_variant_db_contig.h"

namespace kgl = kellerberrin::genome;


// Unconditionally adds a variant to the contig.

bool kgl::ContigDB::addVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_offset_map_.find(variant_ptr->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.

    auto& [contig_offset, offset_ptr] = *result;

    offset_ptr->addVariant(variant_ptr);

  } else {

    // add the new offset.
    std::unique_ptr<OffsetDB> offset_array_ptr(std::make_unique<OffsetDB>());

    offset_array_ptr->addVariant(variant_ptr);

    auto insert_result = contig_offset_map_.try_emplace(variant_ptr->offset(), std::move(offset_array_ptr));

    if (not insert_result.second) {

      ExecEnv::log().error("ContigDB::addVariant(); Could not add variant offset: {} to the genome",
                           variant_ptr->offset());
      return false;

    }

  }

  return true;

}


// Unconditionally adds an offset
bool kgl::ContigDB::addOffset(ContigOffset_t offset, std::unique_ptr<OffsetDB> offset_db) {

  // check variant offsets.
  for (auto const& variant_ptr : offset_db->getVariantArray()) {

    if (offset != variant_ptr->offset()) {

      ExecEnv::log().error( "ContigDB::addOffset; mismatch offset: {}, variant: {}",
                            offset, variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));
      return false;

    }

  }

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_offset_map_.find(offset);

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.

    auto& [contig_offset, offset_ptr] = *result;

    offset_ptr = std::move(offset_db);

  } else {

    auto insert_result = contig_offset_map_.try_emplace(offset, std::move(offset_db));

    if (not insert_result.second) {

      ExecEnv::log().error("ContigDB::addOffset; Could not add variant offset: {} to the genome", offset);
      return false;

    }

  }

  return true;

}

// Unconditionally adds a variant to the contig.
bool kgl::ContigDB::addUniqueUnphasedVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_offset_map_.find(variant_ptr->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.
    const OffsetDBArray& offset_variant_array = result->second->getVariantArray();
    // Check for uniqueness
    for (auto const& offset_variant : offset_variant_array) {

      if (offset_variant->analogous(*variant_ptr)) {

        return true;

      }

    }

    // Variant is unique, so clone, de-phase, and add.
    std::shared_ptr<Variant> unphased_variant = variant_ptr->clonePhase(VariantPhase::UNPHASED);
    result->second->addVariant(unphased_variant);

  } else {
    // add the new offset.
    std::unique_ptr<OffsetDB> offset_array_ptr(std::make_unique<OffsetDB>());
    std::shared_ptr<Variant> unphased_variant = variant_ptr->clonePhase(VariantPhase::UNPHASED);
    offset_array_ptr->addVariant(unphased_variant);
    auto insert_result = contig_offset_map_.try_emplace(variant_ptr->offset(), std::move(offset_array_ptr));

    if (not insert_result.second) {

      ExecEnv::log().error("ContigDB::addUniqueUnphasedVariant(); Could not add variant offset: {} to the genome",
                           variant_ptr->offset());
      return false;

    }

  }

  return true;

}



// Counts the variants in a contug.
size_t kgl::ContigDB::variantCount() const {


  size_t variant_count = 0;

  for (auto const &offset_variant_vector : getMap()) {

    variant_count += offset_variant_vector.second->getVariantArray().size();

  }

  return variant_count;

}


// Creates a copy of the contig that only contains variants passing the filter condition.
std::shared_ptr<kgl::ContigDB> kgl::ContigDB::filterVariants(const VariantFilter &filter) const {

  std::shared_ptr<ContigDB> filtered_contig_ptr(std::make_shared<ContigDB>(contigId()));

  for (auto const&[offset, offset_ptr] : getMap()) {

    std::unique_ptr<OffsetDB> filtered_offset(std::make_unique<OffsetDB>());
    filtered_offset->setVariantArray(offset_ptr->getVariantArray());
    filtered_offset->inSituFilter(filter);
    if (not filtered_offset->getVariantArray().empty()) {

      if (not filtered_contig_ptr->addOffset(offset, std::move(filtered_offset))) {

        ExecEnv::log().error("ContigDB::filterVariants; Problem adding variant at offset: {}, to contig: {}",
                             offset, contigId());

      }

    }

  }

  return filtered_contig_ptr;

}



// Filters variants inSitu
// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
std::pair<size_t, size_t> kgl::ContigDB::inSituFilter(const VariantFilter &filter) {

  std::pair<size_t, size_t> contig_count{0, 0};
  for (auto& [offset, offset_ptr] : contig_offset_map_) {

    auto const offset_count = offset_ptr->inSituFilter(filter);
    contig_count.first += offset_count.first;
    contig_count.second += offset_count.second;

  } // for all offsets

  // Finally, remove all empty offsets.
  auto it = contig_offset_map_.begin();
  while(it != contig_offset_map_.end()) {

    if (it->second->getVariantArray().empty()) {

      it = contig_offset_map_.erase(it);

    } else {

      ++it;

    }

  }

  return contig_count;

}


std::optional<kgl::OffsetDBArray> kgl::ContigDB::findOffsetArray(ContigOffset_t offset) const {

  auto result = contig_offset_map_.find(offset);

  if (result != contig_offset_map_.end()) {

    OffsetDBArray variant_array = result->second->getVariantArray();

    return variant_array;

  } else {

    return std::nullopt;

  }

}



std::optional<std::shared_ptr<const kgl::Variant>> kgl::ContigDB::findVariant(const Variant& variant) const {

  auto result = contig_offset_map_.find(variant.offset());

  if (result != contig_offset_map_.end()) {

    auto& [offset, offset_ptr] = *result;

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (variant.analogous(*variant_ptr)) {

        return variant_ptr;

      }

    }

    return std::nullopt;


  } else {

    return std::nullopt;

  }

}

// This search algorithm has n^2 complexity.
// The variants in the template contig must be unique regardless of phase.
std::shared_ptr<kgl::ContigDB> kgl::ContigDB::findContig(const std::shared_ptr<const ContigDB>& template_contig) const {

  std::shared_ptr<ContigDB> found_contig_ptr(std::make_shared<ContigDB>(template_contig->contigId()));

  for (auto const& [offset, offset_ptr] : template_contig->getMap()) {

    auto result = contig_offset_map_.find(offset);
    if (result != contig_offset_map_.end()) {

      // Create a set of allele hashs to search.
      std::unordered_set<std::string> search_hash;
      for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

        search_hash.insert(variant_ptr->alleleHash());

      }

      // Search the set of hashs.
      auto const& [this_offset, this_offset_ptr] = *result;
      for (auto const& this_variant_ptr : this_offset_ptr->getVariantArray()) {

        auto result = search_hash.find(this_variant_ptr->alleleHash());
        if (result != search_hash.end()) {

          if (not found_contig_ptr->addVariant(this_variant_ptr)) {

            ExecEnv::log().error( "ContigDB::findContig; cannot add variant: {}",
                                  this_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

          }

        }

      }

    } // if this offset

  } // for all template offset

  return found_contig_ptr;

}


std::pair<size_t, size_t> kgl::ContigDB::validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const {

  std::pair<size_t, size_t> contig_count{0, 0};

  std::shared_ptr<const DNA5SequenceContig> contig_sequence_ptr = contig_db_ptr->sequence_ptr();

  for (auto const&[offset, variant_vector] : getMap()) {

    contig_count.first += variant_vector->getVariantArray().size();

    if (offset >= contig_sequence_ptr->length()) {

      ExecEnv::log().error("ContigDB::validate,  Variant offset: {} exceeds total contig: {} size: {}", offset,
                           contig_db_ptr->contigId(), contig_sequence_ptr->length());
      continue;

    }

    for (auto const &variant_ptr : variant_vector->getVariantArray()) {

      if (not variant_ptr) {

        ExecEnv::log().error("ContigDB::validate, Unknown variant: {}",
                             variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));
        continue;

      }

      if (contig_sequence_ptr->subSequence(variant_ptr->referenceOffset(), variant_ptr->reference().length()) == variant_ptr->reference()) {

        ++contig_count.second;

      } else {

        ExecEnv::log().error(" ContigDB::validate, Mismatch, at Contig Offset: {} Sequence is: {}, Variant Reference Sequence is: {}",
                               variant_ptr->offset(),
                               contig_sequence_ptr->subSequence(variant_ptr->offset(),variant_ptr->reference().length()).getSequenceAsString(),
                               variant_ptr->reference().getSequenceAsString());

      }

    }

  }

  return contig_count;

}

bool kgl::ContigDB::getSortedVariants(VariantPhase phase,
                                      ContigOffset_t start,
                                      ContigOffset_t end,
                                      OffsetVariantMap& variant_map) const {


  auto lower_bound = contig_offset_map_.lower_bound(start);
  auto upper_bound = contig_offset_map_.upper_bound(end-1); //  [start, end)

  // If there is a prior variant that overlaps the start address, then push this onto the variant map.

  if (lower_bound != contig_offset_map_.end() and lower_bound != contig_offset_map_.begin()) {

    auto previous_offset_ptr = std::prev(lower_bound);

    const OffsetDBArray& previous_offset_variants = previous_offset_ptr->second->getVariantArray();

    for (const auto& variant : previous_offset_variants) {

      if (variant->phaseId() == phase or phase == VariantPhase::UNPHASED) {

        if (variant->offset() + variant->referenceSize() > start) {

          variant_map.emplace(variant->offset(), variant);

        }

      }

    }

  }



  for (auto it = lower_bound; it != upper_bound; ++it) {


    const OffsetDBArray& previous_offset_variants = it->second->getVariantArray();

    for (const auto& variant : previous_offset_variants) {

      if (variant->phaseId() == phase or phase == VariantPhase::UNPHASED) {

        variant_map.emplace(it->first, variant);

      }

    }

  }

  checkUpstreamDeletion(variant_map);

  return true;

}

void kgl::ContigDB::checkUpstreamDeletion(OffsetVariantMap& variant_map) const {

  for (auto iter = variant_map.begin(); iter != variant_map.end(); ++iter) {

    if (iter == variant_map.begin()) continue;

    int delete_size = std::prev(iter)->second->referenceSize() - std::prev(iter)->second->alternateSize();

    int offset_gap = iter->second->offset() - std::prev(iter)->second->offset();

    if (delete_size >= offset_gap) {

//      ExecEnv::log().info("ContigDB::checkUpstreamDeletion(), Upstream deletion detected: {}",
//                           std::prev(iter)->second->output(' ', VariantOutputIndex::START_0_BASED, true));

//      ExecEnv::log().info("ContigDB::checkUpstreamDeletion(), Downstream variant removed from mutation: {}",
//                           iter->second->output(' ', VariantOutputIndex::START_0_BASED, true));

      iter = variant_map.erase(iter);

      iter = std::prev(iter); // reset back to the previous iterator.

    }

  }

}


std::shared_ptr<kgl::ContigDB> kgl::ContigDB::subset(ContigOffset_t start, ContigOffset_t end) const {

  std::shared_ptr<ContigDB> subset_contig(std::make_shared<ContigDB>(contigId()));

  auto lower_bound = contig_offset_map_.lower_bound(start);
  auto upper_bound = contig_offset_map_.upper_bound(end-1); //  [start, end)

  for (auto it = lower_bound; it != upper_bound; ++it) {

    auto const& [offset, offset_ptr] = *it;

    auto const& variant_array = offset_ptr->getVariantArray();

    for (auto const& variant_ptr : variant_array) {

      if (not subset_contig->addVariant(variant_ptr)) {

        ExecEnv::log().error("ContigDB::subset; problem inserting variant into subset [{}, {}) of contig: {}", start, end, contigId());

      }

    }

  }

  return subset_contig;

}


std::shared_ptr<kgl::ContigDB> kgl::ContigDB::intersection(const std::shared_ptr<const ContigDB>& contig, VariantEquality variant_equality) const {

  std::shared_ptr<ContigDB> intersection_contig(std::make_shared<ContigDB>(contigId()));

  for (auto const& [offset, offset_ptr] : getMap()) {

    auto offset_opt = contig->findOffsetArray(offset);

    if (offset_opt) {

      // Generate a set of unique variants at this offset.
      std::set<std::string> unique_hash_set;
      for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

        switch(variant_equality) {

          case VariantEquality::PHASED:
            unique_hash_set.insert(variant_ptr->allelePhaseHash());
            break;

          case VariantEquality::UNPHASED:
            unique_hash_set.insert(variant_ptr->alleleHash());
            break;

          default:
            break;

        }

      }

      for (auto const& variant_ptr : offset_opt.value()) {

        switch(variant_equality) {

          case VariantEquality::PHASED:
            if (unique_hash_set.find(variant_ptr->allelePhaseHash()) != unique_hash_set.end()) {

              if (not intersection_contig->addVariant(variant_ptr)) {

                ExecEnv::log().error( "ContigDB::intersection, error could not insert variant: {}",
                                      variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

              }

            }
            break;

          case VariantEquality::UNPHASED:
            if (unique_hash_set.find(variant_ptr->alleleHash()) != unique_hash_set.end()) {

              if (not intersection_contig->addVariant(variant_ptr)) {

                ExecEnv::log().error( "ContigDB::intersection, error could not insert variant: {}",
                                      variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

              }

            }
            break;

          default:
            break;

        }

      }

    }

  }

  return intersection_contig;

}


std::shared_ptr<kgl::ContigDB> kgl::ContigDB::complement(const std::shared_ptr<const ContigDB>& contig, VariantEquality variant_equality) const {

  std::shared_ptr<ContigDB> complement(std::make_shared<ContigDB>(contigId()));

  for (auto const& [offset, offset_ptr] : getMap()) {

    auto offset_opt = contig->findOffsetArray(offset);

    if (offset_opt) {

      // Generate a set of unique variants at this offset.
      std::set<std::string> unique_hash_set;
      for (auto const& variant_ptr : offset_opt.value()) {

        switch (variant_equality) {

          case VariantEquality::PHASED:
            unique_hash_set.insert(variant_ptr->allelePhaseHash());
            break;

          case VariantEquality::UNPHASED:
            unique_hash_set.insert(variant_ptr->alleleHash());
            break;

          default:
            break;

        }

      }

      for (auto const &variant_ptr : offset_ptr->getVariantArray()) {

        switch(variant_equality) {

          case VariantEquality::PHASED:
            if (unique_hash_set.find(variant_ptr->allelePhaseHash()) == unique_hash_set.end()) {

              if (not complement->addVariant(variant_ptr)) {

                ExecEnv::log().error( "ContigDB::intersection, error could not insert variant: {}",
                                      variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

              }

            }
            break;

          case VariantEquality::UNPHASED:
            if (unique_hash_set.find(variant_ptr->alleleHash()) == unique_hash_set.end()) {

              if (not complement->addVariant(variant_ptr)) {

                ExecEnv::log().error( "ContigDB::intersection, error could not insert variant: {}",
                                      variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

              }

            }
            break;

          default:
            break;

        }

      }

    } else { //add all variants if offset not found.

      for (auto const &variant_ptr : offset_ptr->getVariantArray()) {

        if (not complement->addVariant(variant_ptr)) {

          ExecEnv::log().error( "ContigDB::intersection, error could not insert variant: {}",
                                variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        }

      }

    }

  }

  return complement;

}
