//
// Created by kellerberrin on 22/04/18.
//


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// HomologousVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::HomologousVariant>
kgl::HomologousVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::HomologousVariant> filtered_contig_ptr = deepCopy();
  // Complements the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
  auto predicate = [&](const OffsetVariantMap::const_iterator& it) { return not it->second->filterVariant(filter); };
  predicateIterableDelete(filtered_contig_ptr->offset_variant_map_,  predicate);

  return filtered_contig_ptr;

}

// This function will insert multiple different variants for contig offset in a std::multimap
bool kgl::HomologousVariant::addVariant(std::shared_ptr<const Variant>& variant_ptr) {

  auto result = offset_variant_map_.equal_range(variant_ptr->offset());

  for (auto it = result.first; it != result.second; ++it) {

    if (variant_ptr->equivalent(*it->second)) {

      std::string variant_str = variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false);
      ExecEnv::log().warn("addVariant() fails, Equivalent variant already exists:\n{}", variant_str);
      return false;

    }

  }

  offset_variant_map_.insert(std::make_pair(variant_ptr->offset(), variant_ptr));

  return true;

}

// Returns true if variant found and erased.
bool kgl::HomologousVariant::eraseVariant(std::shared_ptr<const Variant>& variant_ptr) {

  auto result = offset_variant_map_.equal_range(variant_ptr->offset());

  for (auto it = result.first; it != result.second; ++it) {

    if (variant_ptr->equivalent(*it->second)) {

      // C++ 11 erase() returns next iterator to next valid (or end())
      offset_variant_map_.erase(it);
      return true;

    }

  }

  return false;

}


// Always use deep copy when modifying this object.
std::shared_ptr<kgl::HomologousVariant> kgl::HomologousVariant::deepCopy() const {

  std::shared_ptr<HomologousVariant> copy(std::make_shared<HomologousVariant>(phaseId()));

  for (auto variant : getMap()) {

    copy->addVariant(variant.second);

  }

  return copy;

}


bool kgl::HomologousVariant::getSortedVariants(ContigOffset_t start, ContigOffset_t end, OffsetVariantMap& variant_map) const {

  auto lower_bound = offset_variant_map_.lower_bound(start);
  auto upper_bound = offset_variant_map_.upper_bound(end-1); //  [start, end)

  // If there is a prior variant that overlaps the start address, then push this onto the variant map.
  if (lower_bound != offset_variant_map_.end() and lower_bound != offset_variant_map_.begin()) {

    auto previous_variant_ptr = lower_bound;

    --previous_variant_ptr;

    if (previous_variant_ptr->second->offset() + previous_variant_ptr->second->reference_size() > start) {

      variant_map.insert(std::pair<ContigOffset_t , std::shared_ptr<const Variant>>(previous_variant_ptr->first, previous_variant_ptr->second));

    }

  }

  for (auto it = lower_bound; it != upper_bound; ++it) {

    variant_map.insert(std::pair<ContigOffset_t , std::shared_ptr<const Variant>>(it->first, it->second));

  }

  return true;

}
