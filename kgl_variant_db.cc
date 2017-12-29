//
// Created by kellerberrin on 31/10/17.
//

#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::ContigVariant> filtered_contig_ptr = deepCopy();
  // Complements the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
  auto predicate = [&](const OffsetVariantMap::const_iterator& it) { return not it->second->filterVariant(filter); };
  predicateIterableDelete(filtered_contig_ptr->offset_variant_map_,  predicate);

  return filtered_contig_ptr;

}

// This function will insert multiple variants for contig offset in a std::multimap
void kgl::ContigVariant::addVariant(std::shared_ptr<const Variant>& variant_ptr) {

  offset_variant_map_.insert(std::make_pair(variant_ptr->contigOffset(), variant_ptr));

}

// Always use deep copy when modifying this object.
std::shared_ptr<kgl::ContigVariant> kgl::ContigVariant::deepCopy() const {

  std::shared_ptr<ContigVariant> copy(std::make_shared<ContigVariant>(contigId()));

  for (auto variant : getMap()) {

    copy->addVariant(variant.second);

  }

  return copy;

}


bool kgl::ContigVariant::getSortedVariants(ContigOffset_t start, ContigOffset_t end, OffsetVariantMap& variant_map) const {

  auto lower_bound = offset_variant_map_.lower_bound(start);
  auto upper_bound = offset_variant_map_.upper_bound(end-1); //  [start, end)

  for (auto it = lower_bound; it != upper_bound; ++it) {

    variant_map.insert(std::pair<ContigOffset_t , std::shared_ptr<const Variant>>(it->first, it->second));

  }

  return true;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GenomeVariant::addContigVariant(std::shared_ptr<kgl::ContigVariant>& contig_variant) {

  auto result = genome_variant_map_.insert(std::make_pair(contig_variant->contigId(), contig_variant));

  return result.second;

}

bool kgl::GenomeVariant::getContigVariant(const ContigId_t& contig_id,
                                          std::shared_ptr<ContigVariant>& contig_variant) const {

  auto contig = genome_variant_map_.find(contig_id);

  if (contig != genome_variant_map_.end()) {

    contig_variant = contig->second;
    return true;

  } else {

    contig_variant = nullptr;
    return false;

  }

}


bool kgl::GenomeVariant::addVariant(std::shared_ptr<const Variant> variant) {

  std::shared_ptr<ContigVariant> contig_variant;
  if (not getContigVariant(variant->contigId(), contig_variant)) {

    ExecEnv::log().error("Contig: {} not found, variant: {}",
                         variant->contigId(), variant->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  contig_variant->addVariant(variant);

  return true;

}


std::shared_ptr<kgl::GenomeVariant> kgl::GenomeVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::GenomeVariant> filtered_genome_ptr(std::make_shared<kgl::GenomeVariant>(genomeId()));

  filtered_genome_ptr->attributes().insertAttribute("filter", filter.filterName());

  ExecEnv::log().info("Applying filter: {}", filter.filterName());
  for (const auto& contig_variant : genome_variant_map_) {

    std::shared_ptr<kgl::ContigVariant> filtered_contig = contig_variant.second->filterVariants(filter);
    filtered_genome_ptr->addContigVariant(filtered_contig);
    ExecEnv::log().vinfo("Contig: {} has: {} filtered variants", contig_variant.first, filtered_contig->variantCount());

  }

  return filtered_genome_ptr;

}

// Creates an empty genome variant with the same contig structure as the genome database.
std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::emptyGenomeVariant(const GenomeId_t& genome_id,
                                       const std::shared_ptr<const GenomeDatabase>& genome_db) {


  std::shared_ptr<GenomeVariant> empty_genome_variant(std::make_shared<GenomeVariant>(genome_id));

  for (auto contig_db : genome_db->getMap()) {

    std::shared_ptr<ContigVariant> contig_variant(std::make_shared<ContigVariant>(contig_db.first));
    if (not empty_genome_variant->addContigVariant(contig_variant)) {

      ExecEnv::log().error("emptyGenomeVariant(), could not add contig variant: {}", contig_db.first);

    }

  }

  return empty_genome_variant;

}




size_t kgl::GenomeVariant::size() const {

  size_t total_variants = 0;
  for (auto contig_variant : genome_variant_map_) {

    total_variants += contig_variant.second->size();

  }

  return total_variants;

}



std::string kgl::GenomeVariant::output(char field_delimter, VariantOutputIndex output_index, bool detail) const {

  std::stringstream ss;
  for (const auto& contig_variant : genome_variant_map_) {

    for (const auto& variant : contig_variant.second->getMap()) {

      ss << variant.second->output(field_delimter, output_index, detail);

    }

  }

  return ss.str();

}


void kgl::GenomeVariant::getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  for (const auto& contig_variant : genome_variant_map_) {

    for (const auto& variant : contig_variant.second->getMap()) {

      variant_vector.push_back(variant.second);

    }

  }

}

// Always use deep copy when modifying this object.
std::shared_ptr<kgl::GenomeVariant> kgl::GenomeVariant::deepCopy() const {

  std::shared_ptr<GenomeVariant> copy(std::make_shared<GenomeVariant>(genomeId()));
  copy->attributes(attributes());

  for (auto contig : getMap()) {

    std::shared_ptr<ContigVariant> copy_contig(std::make_shared<ContigVariant>(contig.second->contigId()));
    copy->addContigVariant(copy_contig);

  }

  std::vector<std::shared_ptr<const Variant>> variant_vector;
  getVariants(variant_vector);

  for (auto variant : variant_vector) {

    copy->addVariant(variant);

  }

  return copy;

}

bool kgl::GenomeVariant::getSortedVariants(ContigId_t contig_id,
                                           ContigOffset_t start,
                                           ContigOffset_t end,
                                           OffsetVariantMap& variant_map) const {

  std::shared_ptr<ContigVariant> contig_variant_ptr;
  if (not getContigVariant(contig_id, contig_variant_ptr)) {

    ExecEnv::log().error("Contig Id: {} not found in Genome Variant: {}", contig_id, genomeId());
    return false;

  }

  return contig_variant_ptr->getSortedVariants(start, end, variant_map);

}


bool kgl::GenomeVariant::getCodingSortedVariants(ContigId_t contig_id,
                                                 ContigOffset_t start,
                                                 ContigOffset_t end,
                                                 OffsetVariantMap& variant_map) const {

  OffsetVariantMap all_variant_map;
  if (not getSortedVariants(contig_id, start, end, all_variant_map)) {

    return false;

  }

  for (auto variant : all_variant_map) {

    if (variant.second->type() == VariantSequenceType::CDS_CODING) {

      variant_map.insert(std::pair<ContigOffset_t,std::shared_ptr<const Variant>>(variant.first, variant.second));

    }

  }

  return true;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::PopulationVariant::addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant) {

  auto result = population_variant_map_.insert(std::pair<GenomeId_t, std::shared_ptr<const GenomeVariant>>(genome_variant->genomeId(), genome_variant));

  return result.second;

}

bool kgl::PopulationVariant::getGenomeVariant(const GenomeId_t& genome_id, std::shared_ptr<const GenomeVariant>& genome_variant) {

  auto result = population_variant_map_.find(genome_id);

  if (result != population_variant_map_.end()) {

    genome_variant = result->second;
    return true;

  } else {

    genome_variant = nullptr;
    return false;

  }

}
