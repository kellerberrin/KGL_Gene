//
// Created by kellerberrin on 25/12/20.
//


#include "kgl_variant_db_contig.h"
#include "kgl_variant_filter_db_variant.h"

#include <ranges>

namespace kgl = kellerberrin::genome;

// Use this to copy the object.
std::shared_ptr<kgl::ContigDB> kgl::ContigDB::deepCopy() const {

  // Can use the shallow filter because all variants are copied across.
  return viewFilter(TrueFilter());

}


// Unconditionally adds a variant to the contig_ref_ptr.
bool kgl::ContigDB::addVariant(const std::shared_ptr<const Variant> &variant_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(lock_contig_mutex_);

  auto [insert_iter, insert_result] = contig_offset_map_.try_emplace(variant_ptr->offset());
  if (insert_result) {

    // Add the new offset.
    insert_iter->second = std::make_shared<OffsetDB>();

  }

  insert_iter->second->addVariant(variant_ptr);

  return true;

}


// Unconditionally adds an offset
bool kgl::ContigDB::addOffset(ContigOffset_t offset, std::unique_ptr<OffsetDB> offset_db) {

  // check variant offsets.
  for (auto const& variant_ptr : offset_db->getVariantArray()) {

    if (offset != variant_ptr->offset()) {

      ExecEnv::log().error( "ContigDB::addOffset; mismatch offset: {}, variant: {}", offset, variant_ptr->HGVS());
      return false;

    }

  }

  // Lock this function to concurrent access.
  std::scoped_lock lock(lock_contig_mutex_);

  auto [insert_iter, insert_result] = contig_offset_map_.try_emplace(offset);
  if (insert_result) {

    // Add the new offset.
    insert_iter->second = std::shared_ptr<OffsetDB>(std::move(offset_db));

  } else {

    // Note that this overwrites any offset that already exists at this offset.
    insert_iter->second = std::shared_ptr<OffsetDB>(std::move(offset_db));

  }

  return true;

}


// Counts the variants in a contig.
size_t kgl::ContigDB::variantCount() const {

  return std::ranges::fold_left(getMap() | std::views::values,
                                size_t{0},
                                [](size_t sum, const auto& offset_ptr) { return sum + offset_ptr->getVariantArray().size(); });

}


// Creates a copy of the contig_ref_ptr that only contains variants passing the filter condition.
// Note that we delete any empty offsets.
std::unique_ptr<kgl::ContigDB> kgl::ContigDB::viewFilter(const BaseFilter &filter) const {

  // Only contig_ref_ptr filter is implemented at this level.
  if (filter.filterType() == FilterBaseType::CONTIG_FILTER) {

    // The filter is cloned before use because filters may maintain (mutable) state
    // that must not transfer between filter invocations or between concurrent threads.
    std::shared_ptr<const FilterContigs> contig_filter = std::dynamic_pointer_cast<const FilterContigs>(filter.clone());
    auto filtered_contig_ptr = contig_filter->applyFilter(*this);
    filtered_contig_ptr->trimEmpty();  // Remove any empty offsets.
    return filtered_contig_ptr;

  }

  // All other filters.
  // Filter the offsets.
  std::unique_ptr<ContigDB> filtered_contig_ptr(std::make_unique<ContigDB>(contigId()));
  for (const auto& [offset, offset_ptr] : getMap()) {

    auto filtered_offset_ptr = offset_ptr->viewFilter(filter);
    if (not filtered_contig_ptr->addOffset(offset, std::move(filtered_offset_ptr))) {

      ExecEnv::log().error("ContigDB::viewFilter; Problem adding offset: {}, to contig_ref_ptr: {}", offset, contigId());

    }

  }

  filtered_contig_ptr->trimEmpty();  // Remove any empty offsets.
  return filtered_contig_ptr;

}


// Filters inSitu
// Returns a std::pair with .first the reference number of variants, .second the filtered number of variants.
// Note that we delete any empty offsets.
std::pair<size_t, size_t> kgl::ContigDB::selfFilter(const BaseFilter &filter) {

  // Only contig filter is implemented at this level.
  if (filter.filterType() == FilterBaseType::CONTIG_FILTER) {

    // The filter is cloned before use because filters may maintain (mutable) state
    // that must not transfer between filter invocations or between concurrent threads.
    std::shared_ptr<const FilterContigs> contig_filter = std::dynamic_pointer_cast<const FilterContigs>(filter.clone());

    size_t prior_count = variantCount();

    auto filtered_contig_ptr = contig_filter->applyFilter(*this);
    contig_offset_map_ = std::move(filtered_contig_ptr->contig_offset_map_);
    trimEmpty();  // Remove any empty offsets.

    size_t post_count = variantCount();

    return { prior_count, post_count };

  }

  // All other filters.
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

  return std::erase_if(contig_offset_map_, [](const auto& offset_pair) { return offset_pair.second->getVariantArray().empty(); });

}


std::optional<kgl::OffsetDBArray> kgl::ContigDB::findOffsetArray(ContigOffset_t offset) const {

  // Note that the variant array is copied on return (the public API requires the copy).
  auto result = contig_offset_map_.find(offset);

  if (result != contig_offset_map_.end()) {

    return result->second->getVariantArray();

  }

  return std::nullopt;

}


std::pair<size_t, size_t> kgl::ContigDB::validate(const std::shared_ptr<const ContigReference> &contig_db_ptr) const {

  std::pair<size_t, size_t> contig_count{0, 0};

  const DNA5SequenceLinear& contig_sequence = contig_db_ptr->sequence();

  for (auto const& [offset, variant_vector] : getMap()) {

    contig_count.first += variant_vector->getVariantArray().size();

    if (offset >= contig_sequence.length()) {

      ExecEnv::log().error("Variant offset: {} exceeds total contig_ref_ptr: {} size: {}", offset,
                           contig_db_ptr->contigId(), contig_sequence.length());
      continue;

    }

    for (auto const &variant_ptr : variant_vector->getVariantArray()) {

      if (not variant_ptr) {

        ExecEnv::log().error("ContigDB::validate, unexpected NULL variant found");
        continue;

      }

      OpenRightUnsigned contig_ref_interval(variant_ptr->offset(), variant_ptr->offset()+variant_ptr->referenceSize());
      auto contig_ref_opt = contig_sequence.subSequence(contig_ref_interval);
      if (not contig_ref_opt) {

        ExecEnv::log().error("Unable to extract variant reference from contig: {}, contig_ref_ptr interval: {}, variant: {}",
                             contig_db_ptr->contigId(),
                             contig_sequence.interval().toString(),
                             variant_ptr->HGVS());
        continue;

      }

      if (contig_ref_opt.value() == variant_ptr->reference()) {

        ++contig_count.second;

      } else {

        ExecEnv::log().error("Mismatch, at Contig Offset: {} Contig Sequence is: {}, Variant is: {}",
                              variant_ptr->offset(),
                              contig_ref_opt.value().getStringView(),
                              variant_ptr->HGVS());

      }

    }

  }

  return contig_count;

}


// Create an equivalent contig_ref_ptr that has canonical variants, SNP are represented by '1X', Deletes by '1MnD'
// and Inserts by '1MnI'. The population structure is re-created and is not a shallow copy.
std::unique_ptr<kgl::ContigDB> kgl::ContigDB::canonicalContig() const {

  // Create the new contig and populate with canonical variants.
  auto canonical_contig_ptr = std::make_unique<ContigDB>(contigId());
  processAll([&canonical_contig_ptr](const std::shared_ptr<const Variant>& variant_ptr) -> bool {

    std::shared_ptr<const Variant> canonical_variant_ptr = variant_ptr->cloneCanonical();

    if (not canonical_contig_ptr->addVariant(canonical_variant_ptr)) {

      ExecEnv::log().error("ContigDB::canonicalContig; Could NOT add canonical variant: {}", canonical_variant_ptr->HGVS());

    }

    return true;

  });

  return canonical_contig_ptr;

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