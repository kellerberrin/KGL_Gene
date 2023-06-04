//
// Created by kellerberrin on 25/12/20.
//


#include "kgl_variant_db_contig.h"
#include "kgl_variant_filter.h"

namespace kgl = kellerberrin::genome;

// Use this to copy the object.
std::shared_ptr<kgl::ContigDB> kgl::ContigDB::deepCopy() const {

  // Can use the shallow filter because all variants are copied across.
  return viewFilter(TrueFilter());

}

// Unconditionally adds a variant to the contig.
bool kgl::ContigDB::addVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(lock_contig_mutex_);

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
  std::scoped_lock lock(lock_contig_mutex_);

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


// Counts the variants in a contug.
size_t kgl::ContigDB::variantCount() const {


  size_t variant_count{0};

  for (auto const &[offset, variant_vector_ptr] : getMap()) {

    variant_count += variant_vector_ptr->getVariantArray().size();

  }

  return variant_count;

}


// Creates a copy of the contig that only contains variants passing the filter condition.
// Note that we delete any empty offsets.
std::unique_ptr<kgl::ContigDB> kgl::ContigDB::viewFilter(const BaseFilter &filter) const {

  // Only contig filter is implemented at this level.
  std::shared_ptr<const FilterContigs> contig_filter = std::dynamic_pointer_cast<const FilterContigs>(filter.clone());
  if (contig_filter) {

    auto filtered_contig_ptr = contig_filter->applyFilter(*this);
    filtered_contig_ptr->trimEmpty();  // Remove any empty offsets.
    return filtered_contig_ptr;

  }

  // All other filters.
  // Filter the offsets.
  std::unique_ptr<ContigDB> filtered_contig_ptr(std::make_unique<ContigDB>(contigId()));
  for (const auto& [offset, offset_ptr] : getMap()) {

    auto filtered_offset_ptr = offset_ptr->copyFilter(filter);
    if (not filtered_contig_ptr->addOffset(offset, std::move(filtered_offset_ptr))) {

      ExecEnv::log().error("ContigDB::filter; Problem adding offset: {}, to contig: {}", offset, contigId());

    }

  }

  filtered_contig_ptr->trimEmpty();  // Remove any empty offsets.
  return filtered_contig_ptr;

}

// Filters inSitu
// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
// Note that we delete any empty offsets.
std::pair<size_t, size_t> kgl::ContigDB::selfFilter(const BaseFilter &filter) {

  // Only genome filter is implemented at this level.
  std::shared_ptr<const FilterContigs> contig_filter = std::dynamic_pointer_cast<const FilterContigs>(filter.clone());
  if (contig_filter) {

    size_t prior_count = variantCount();

    auto filtered_contig_ptr = contig_filter->applyFilter(*this);
    contig_offset_map_ = std::move(filtered_contig_ptr->contig_offset_map_);
    trimEmpty();  // Remove any empty offsets.

    size_t post_count = variantCount();

    return { prior_count, post_count };

  }

  std::pair<size_t, size_t> contig_count{0, 0};
  for (auto& [offset, offset_ptr] : contig_offset_map_) {

    auto const offset_count = offset_ptr->selfFilter(filter);
    contig_count.first += offset_count.first;
    contig_count.second += offset_count.second;

  } // for all offsets

  trimEmpty();  // Remove any empty offsets.
  return contig_count;

}

// Deletes any empty Offsets, returns number deleted.
size_t kgl::ContigDB::trimEmpty() {

  size_t delete_count{0};
  // Delete empty genomes.
  auto it = contig_offset_map_.begin();
  while (it != contig_offset_map_.end()) {

    if (it->second->getVariantArray().empty()) {

      it = contig_offset_map_.erase(it);
      ++delete_count;

    } else {

      ++it;

    }

  }

  return delete_count;

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
// The variants in the template contig are unique. Variant phase is disregarded.
std::shared_ptr<kgl::ContigDB> kgl::ContigDB::findContig(const std::shared_ptr<const ContigDB>& template_contig) const {

  std::shared_ptr<ContigDB> found_contig_ptr(std::make_shared<ContigDB>(template_contig->contigId()));

  for (auto const& [offset, offset_ptr] : template_contig->getMap()) {

    auto result = contig_offset_map_.find(offset);
    if (result != contig_offset_map_.end()) {

      // Create a set of allele hashs to search.
      std::unordered_set<std::string> search_hash;
      for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

        search_hash.insert(variant_ptr->HGVS());

      }

      // Search the set of hashs.
      auto const& [this_offset, this_offset_ptr] = *result;
      for (auto const& this_variant_ptr : this_offset_ptr->getVariantArray()) {

        auto result = search_hash.find(this_variant_ptr->HGVS());
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

  for (auto const& [offset, variant_vector] : getMap()) {

    contig_count.first += variant_vector->getVariantArray().size();

    if (offset >= contig_sequence_ptr->length()) {

      ExecEnv::log().error("ContigDB::validate,  Variant offset: {} exceeds total contig: {} size: {}", offset,
                           contig_db_ptr->contigId(), contig_sequence_ptr->length());
      continue;

    }

    for (auto const &variant_ptr : variant_vector->getVariantArray()) {

      if (not variant_ptr) {

        ExecEnv::log().error("ContigDB::validate, unexpected NULL variant found");
        continue;

      }

      if (contig_sequence_ptr->subSequence(variant_ptr->offset(), variant_ptr->reference().length()) == variant_ptr->reference()) {

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

bool kgl::ContigDB::processAll(const VariantProcessFunc& objFunc)  const {

  for (auto const& [offset, offset_ptr] : getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (not objFunc(variant_ptr)) {

        ExecEnv::log().error("ContigDB::processAll; Problem executing general purpose function at offset: {}", offset);
        return false;

      }

    }

  }

  return true;

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

    auto const& [offset, offset_db_ptr] = *previous_offset_ptr;

    const OffsetDBArray& previous_offset_variants = offset_db_ptr->getVariantArray();

    for (const auto& variant_ptr : previous_offset_variants) {

      if (variant_ptr->phaseId() == phase or phase == VariantPhase::UNPHASED) {

        if (variant_ptr->offset() + variant_ptr->referenceSize() > start) {

          variant_map.emplace(variant_ptr->offset(), variant_ptr);

        }

      }

    }

  }

  // Get all variants between the lower and upper bounds.
  for (auto it = lower_bound; it != upper_bound; ++it) {

    auto const& [offset, offset_db_ptr] = *it;

    const OffsetDBArray& previous_offset_variants = offset_db_ptr->getVariantArray();

    for (const auto& variant_ptr : previous_offset_variants) {

      if (variant_ptr->phaseId() == phase or phase == VariantPhase::UNPHASED) {

        variant_map.emplace(offset, variant_ptr);

      }

    }

  }

  checkUpstreamDeletion(variant_map);

  return true;

}

void kgl::ContigDB::checkUpstreamDeletion(OffsetVariantMap& variant_map) const {

  for (auto iter = variant_map.begin(); iter != variant_map.end(); ++iter) {

    // Ignore the first entry.
    if (iter == variant_map.begin()) {

      continue;

    }

    auto const& [previous_offset, previous_variant_ptr] = *std::prev(iter);
    auto const& [offset, variant_ptr] = *iter;

    int64_t delete_size = previous_variant_ptr->referenceSize() - previous_variant_ptr->alternateSize();
    int64_t offset_gap = variant_ptr->offset() - previous_variant_ptr->offset();

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


// Find the [lower, upper] offsets of a contig, {0,0} is returned if empty.
std::pair<kgl::ContigOffset_t, kgl::ContigOffset_t> kgl::ContigDB::offsetBounds() const {

  if (contig_offset_map_.empty()) {

    return {0, 0};

  }

  auto const& [lower_bound, lower_offset_db] = *contig_offset_map_.begin();
  auto const& [upper_bound, upper_offset_db] = *contig_offset_map_.rbegin();

  return {lower_bound, upper_bound};

}


std::shared_ptr<kgl::ContigDB> kgl::ContigDB::subset(ContigOffset_t start, ContigOffset_t end) const {

  std::shared_ptr<ContigDB> subset_contig(std::make_shared<ContigDB>(contigId()));

  auto lower_bound = contig_offset_map_.lower_bound(start);
  auto upper_bound = contig_offset_map_.upper_bound(end); //  [start, end)

  while (lower_bound != upper_bound) {

    auto const& [offset, offset_ptr] = *lower_bound;

    auto const& variant_array = offset_ptr->getVariantArray();

    for (auto const& variant_ptr : variant_array) {

      if (not subset_contig->addVariant(variant_ptr)) {

        ExecEnv::log().error("ContigDB::subset; problem inserting variant into subset [{}, {}) of contig: {}", start, end, contigId());

      }

    }

    ++lower_bound;

  }

  return subset_contig;

}


std::unique_ptr<kgl::ContigDB> kgl::ContigDB::setIntersection(const ContigDB& contig_B, VariantEquality variant_equality) const {

  std::unique_ptr<ContigDB> intersection_contig(std::make_unique<ContigDB>(contigId()));

  for (auto const& [offset, offset_ptr] : getMap()) {

    auto result = contig_B.getMap().find(offset);
    if (result != contig_B.getMap().end()) {

      auto const& [offset, offset_db_ptr] = *result;

      auto intersection_offset = offset_ptr->setIntersection(*offset_db_ptr, variant_equality);
      // If not empty add the intersection.
      if (not intersection_offset->getVariantArray().empty()) {

        ContigOffset_t offset = intersection_offset->getVariantArray().front()->offset();
        if (not intersection_contig->addOffset(offset, std::move(intersection_offset))) {

          ExecEnv::log().error("ContigDB::setIntersection; Cannot add offset: {}", offset);

        }

      }

    }

  }

  return intersection_contig;

}


std::unique_ptr<kgl::ContigDB> kgl::ContigDB::setComplement(const ContigDB& contig_B, VariantEquality variant_equality) const {

  std::unique_ptr<ContigDB> complement(std::make_unique<ContigDB>(contigId()));

  for (auto const& [offset, offset_ptr] : getMap()) {

    auto result = contig_B.getMap().find(offset);
    if (result == contig_B.getMap().end()) {

      // Create a copy of this offset by taking the complement of the empty offset.
      auto offset_copy = offset_ptr->setComplement(OffsetDB(), variant_equality);
      if (not complement->addOffset(offset, std::move(offset_copy))) {

        ExecEnv::log().error("ContigDB::setComplement; Cannot add offset: {}", offset);

      }

    } else {

      auto complement_offset = offset_ptr->setComplement(*(result->second), variant_equality);
      // If not empty add the complement.
      if (not complement_offset->getVariantArray().empty()) {

        if (not complement->addOffset(offset, std::move(complement_offset))) {

          ExecEnv::log().error("ContigDB::setComplement; Cannot add offset: {}", offset);

        }

      }

    }

  }

  return complement;

}


std::unique_ptr<kgl::ContigDB> kgl::ContigDB::setUnion(const ContigDB& contig_B, VariantEquality variant_equality) const {

  std::unique_ptr<ContigDB> union_contig(std::make_unique<ContigDB>(contigId()));

  for (auto const& [offset, offset_ptr] : getMap()) {

    auto result = contig_B.getMap().find(offset);
    if (result == contig_B.getMap().end()) {

      // Create a unique copy of this offset by taking the union_contig of the empty offset.
      auto offset_copy = offset_ptr->setUnion(OffsetDB(), variant_equality);
      if (not union_contig->addOffset(offset, std::move(offset_copy))) {

        ExecEnv::log().error("ContigDB::setUnion; Cannot add offset: {}", offset);

      }

    } else {

      auto union_offset = offset_ptr->setUnion(*(result->second), variant_equality);
      // If not empty add the union_contig.
      if (not union_offset->getVariantArray().empty()) {

        if (not union_contig->addOffset(offset, std::move(union_offset))) {

          ExecEnv::log().error("ContigDB::setUnion; Cannot add offset: {}", offset);

        }

      }

    }

  }

  return union_contig;

}
