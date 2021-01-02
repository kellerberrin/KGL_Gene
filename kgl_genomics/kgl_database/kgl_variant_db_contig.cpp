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


// Unconditionally adds a variant to the contig.
bool kgl::ContigDB::addUniqueUnphasedVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = contig_offset_map_.find(variant_ptr->offset());

  if (result != contig_offset_map_.end()) {
    // Variant offset exists.
    OffsetDBArray offset_variant_array = result->second->getVariantArray();
    // Check for uniqueness
    for (auto const& offset_variant : offset_variant_array) {

      if (offset_variant->analogous(*variant_ptr)) {

        return true;

      }

    }

    // Variant is unique, so clone, de-phase, and add.
    std::shared_ptr<Variant> unphased_variant = variant_ptr->clone();
    unphased_variant->updatePhaseId(VariantSequence::UNPHASED);
    result->second->addVariant(unphased_variant);

  } else {
    // add the new offset.
    std::unique_ptr<OffsetDB> offset_array_ptr(std::make_unique<OffsetDB>());
    std::shared_ptr<Variant> unphased_variant = variant_ptr->clone();
    unphased_variant->updatePhaseId(VariantSequence::UNPHASED);
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

  for (auto const&[offset, variant_vector] : getMap()) {

    for (auto const &variant_ptr : variant_vector->getVariantArray()) {

      if (filter.applyFilter(*variant_ptr)) {

        if (not filtered_contig_ptr->addVariant(variant_ptr)) {

          ExecEnv::log().error("ContigDB::filterVariants; Problem adding variant at offset: {}, to contig: {}",
                               offset, contigId());

        }

      }

    }

  }

  return filtered_contig_ptr;

}



// Filters variants inSitu
// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
std::pair<size_t, size_t> kgl::ContigDB::inSituFilter(const VariantFilter &filter) {

  std::pair<size_t, size_t> contig_count{0, 0};
  for (auto& [offset, variant_vector] : contig_offset_map_) {

    OffsetDBArray filtered_variants;
    const OffsetDBArray& variant_array = variant_vector->getVariantArray();
    contig_count.first += variant_array.size();

    for (auto const &variant_ptr : variant_array) {

      if (filter.applyFilter(*variant_ptr)) {

        filtered_variants.push_back(variant_ptr);
        ++contig_count.second;

      }

      filtered_variants.push_back(variant_ptr);

    }

    // Replace with filtered variants at this offset.
    variant_vector->replaceVector(std::move(filtered_variants));

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

    OffsetDBArray variant_array = result->second->getVariantArray();

    for (auto const& variant_ptr : variant_array) {

      if (variant.analogous(*variant_ptr)) {

        return variant_ptr;

      }

    }

    return std::nullopt;


  } else {

    return std::nullopt;

  }

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

bool kgl::ContigDB::getSortedVariants(PhaseId_t phase,
                                      ContigOffset_t start,
                                      ContigOffset_t end,
                                      OffsetVariantMap& variant_map) const {


  auto lower_bound = contig_offset_map_.lower_bound(start);
  auto upper_bound = contig_offset_map_.upper_bound(end-1); //  [start, end)

  // If there is a prior variant that overlaps the start address, then push this onto the variant map.

  if (lower_bound != contig_offset_map_.end() and lower_bound != contig_offset_map_.begin()) {

    auto previous_offset_ptr = std::prev(lower_bound);

    OffsetDBArray previous_offset_variants = previous_offset_ptr->second->getVariantArray();

    for (const auto& variant : previous_offset_variants) {

      if (variant->phaseId() == phase) {

        if (variant->offset() + variant->referenceSize() > start) {

          variant_map.emplace(variant->offset(), variant);

        }

      }

    }

  }



  for (auto it = lower_bound; it != upper_bound; ++it) {


    OffsetDBArray previous_offset_variants = it->second->getVariantArray();

    for (const auto& variant : previous_offset_variants) {

      if (variant->phaseId() == phase) {

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



