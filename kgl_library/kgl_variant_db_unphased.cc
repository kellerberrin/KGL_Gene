//
// Created by kellerberrin on 23/04/18.
//


#include "kgl_variant_db_unphased.h"
#include "kgl_patterns.h"


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a contig.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::UnphasedContig::addVariant(std::shared_ptr<const Variant> variant) {

  auto result = contig_offset_map_.find(variant->offset());

  if (result != contig_offset_map_.end()) {

    result->second.push_back(variant);

    return true;

  } else {

    std::pair<ContigOffset_t, std::vector<std::shared_ptr<const Variant>>> new_offset;
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
bool kgl::UnphasedContig::isHomozygous(const std::vector<std::shared_ptr<const Variant>>& variant_vector) {

  bool is_homozygous = true;
  for (auto variant : variant_vector) {

    if (not variant_vector.front()->equivalent(*variant)) {

      is_homozygous = false;
      break;

    }

  }

  return is_homozygous;

}


std::shared_ptr<kgl::UnphasedContig> kgl::UnphasedContig::removeConflictingVariants() const {

  std::shared_ptr<UnphasedContig> filtered_contig_ptr(std::make_shared<UnphasedContig>(contigId()));

  for (auto& variant_vector : getMap()) {

    if (not variant_vector.second.empty()) {

      filtered_contig_ptr->addVariant(variant_vector.second.front());

    } else {

      ExecEnv::log().warn("removeConflictingVariants(); empty variant vector in unphased variant contig: {}, offset: {}",
                          contig_id_, variant_vector.first);
    }

  }

  return filtered_contig_ptr;

}


std::shared_ptr<kgl::UnphasedContig> kgl::UnphasedContig::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::UnphasedContig> filtered_contig_ptr(std::make_shared<kgl::UnphasedContig>(contigId()));

  // Complements the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
  auto predicate = [&](const UnphasedVariantVector::const_iterator& it) { return not (*it)->filterVariant(filter); };

  for (auto offset_vector : getMap()) {

    std::vector<std::shared_ptr<const Variant>> copy_offset_vector = offset_vector.second;

    predicateIterableDelete(copy_offset_vector,  predicate);

    for (auto variant : copy_offset_vector) {

      filtered_contig_ptr->addVariant(variant);

    }

  }

  return filtered_contig_ptr;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a genome.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::UnphasedGenome::addVariant(std::shared_ptr<const Variant> variant) {

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


bool kgl::UnphasedGenome::addContig(std::shared_ptr<UnphasedContig> contig_ptr) {

  std::pair<ContigId_t, std::shared_ptr<UnphasedContig>> add_contig(contig_ptr->contigId(), contig_ptr);
  auto result = contig_map_.insert(add_contig);

  if (not result.second) {

    ExecEnv::log().error("UnphasedGenome::addContig(); could not add contig: {} to the genome", contig_ptr->contigId());

  }

  return result.second;

}


std::shared_ptr<kgl::UnphasedGenome> kgl::UnphasedGenome::removeConflictingVariants() const {

  std::shared_ptr<UnphasedGenome> filtered_genome_ptr(std::make_shared<UnphasedGenome>(genomeId()));

  for (auto contig : getMap()) {

    filtered_genome_ptr->addContig(contig.second->removeConflictingVariants());

  }

  return filtered_genome_ptr;

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
    filtered_genome_ptr->addContig(filtered_contig);
    ExecEnv::log().vinfo("Contig: {} has: {} filtered variants", contig_variant.first, filtered_contig->variantCount());

  }

  return filtered_genome_ptr;

}

