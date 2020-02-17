//
// Created by kellerberrin on 23/04/18.
//


#include <memory>
#include "kel_patterns.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const kgl::PhaseId_t kgl::ContigVariant::HAPLOID_HOMOLOGOUS_INDEX;

std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<ContigVariant> filtered_contig_ptr(std::make_shared<ContigVariant>(contigId(), ploidy()));

  size_t index = 0;
  for (auto homologous : ploidy_vector_) {

    filtered_contig_ptr->ploidy_vector_[index] = homologous->filterVariants(filter);
    index++;

  }

  return filtered_contig_ptr;

}

// This function will insert multiple different variants for contig offset in a std::multimap
bool kgl::ContigVariant::addVariant(std::shared_ptr<const Variant>& variant_ptr) {


  if (variant_ptr->phaseId() >= ploidy()) {

    ExecEnv::log().error("ContigVariant::addVariant(); Contig ploidy: {}, Variant has incompatible phasing: {}",
                         ploidy(), variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  return ploidy_vector_[variant_ptr->phaseId()]->addVariant(variant_ptr);

}

// Returns true if variant found and erased.
bool kgl::ContigVariant::eraseVariant(std::shared_ptr<const Variant>& variant_ptr) {

  if (variant_ptr->phaseId() >= ploidy()) {

    ExecEnv::log().error("ContigVariant::eraseVariant(); Contig ploidy: {}, Variant has incompatible phasing: {}",
                         ploidy(), variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  return ploidy_vector_[variant_ptr->phaseId()]->eraseVariant(variant_ptr);

}


// Always use deep copy when modifying this object.
std::shared_ptr<kgl::ContigVariant> kgl::ContigVariant::deepCopy() const {

  std::shared_ptr<ContigVariant> copy(std::make_shared<ContigVariant>(contigId(), ploidy()));

  for (auto homologous : getVector()) {

    for (auto variant : homologous->getMap()) {

      copy->addVariant(variant.second);

    }

  }

  return copy;

}


bool kgl::ContigVariant::getSortedVariants(PhaseId_t phase,
                                           ContigOffset_t start,
                                           ContigOffset_t end,
                                           OffsetVariantMap& variant_map) const {

  if (phase >= ploidy()) {

    ExecEnv::log().error("ContigVariant::getSortedVariants(); Incompatible phase: {}, Contig ploidy: {}", phase, ploidy());
    variant_map.clear();
    return false;

  }

  return ploidy_vector_[phase]->getSortedVariants(start, end, variant_map);

}


size_t kgl::ContigVariant::variantCount() const {

  size_t variant_count = 0;
  for (auto homologous : getVector()) {

    variant_count += homologous->variantCount();

  }

  return variant_count;

}


