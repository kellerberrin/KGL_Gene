//
// Created by kellerberrin on 31/10/17.
//

#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Ploidy constants.
const kgl::PhaseId_t kgl::GenomeVariant::HAPLOID_GENOME;
const kgl::PhaseId_t kgl::GenomeVariant::DIPLOID_GENOME;


bool kgl::GenomeVariant::addContigVariant(std::shared_ptr<kgl::ContigVariant>& contig_variant) {

  if (contig_variant->ploidy() != ploidy()) {

    ExecEnv::log().error("GenomeVariant::addContigVariant(); Ploidy mismatch - genome ploidy: {}, contig ploidy: {}"
                         , ploidy(), contig_variant->ploidy());
    return false;

  }

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

  return contig_variant->addVariant(variant);

}


std::shared_ptr<kgl::GenomeVariant> kgl::GenomeVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::GenomeVariant> filtered_genome_ptr(std::make_shared<kgl::GenomeVariant>(genomeId(), ploidy()));

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
                                       PhaseId_t ploidy,
                                       const std::shared_ptr<const GenomeDatabase>& genome_db) {


  std::shared_ptr<GenomeVariant> empty_genome_variant(std::make_shared<GenomeVariant>(genome_id, ploidy));

  for (auto contig_db : genome_db->getMap()) {

    std::shared_ptr<ContigVariant> contig_variant(std::make_shared<ContigVariant>(contig_db.first, ploidy));
    if (not empty_genome_variant->addContigVariant(contig_variant)) {

      ExecEnv::log().error("emptyGenomeVariant(), could not add contig variant: {}", contig_db.first);

    }

  }

  return empty_genome_variant;

}




size_t kgl::GenomeVariant::variantCount() const {

  size_t total_variants = 0;
  for (auto contig_variant : genome_variant_map_) {

    total_variants += contig_variant.second->variantCount();

  }

  return total_variants;

}



std::string kgl::GenomeVariant::output(char field_delimter, VariantOutputIndex output_index, bool detail) const {

  std::stringstream ss;
  for (const auto& contig_variant : genome_variant_map_) {

    for (const auto& homologous : contig_variant.second->getVector()) {

      for (const auto variant : homologous->getMap()) {

        ss << variant.second->output(field_delimter, output_index, detail);

      }

    }

  }

  return ss.str();

}


void kgl::GenomeVariant::getVariants(std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  for (const auto& contig_variant : genome_variant_map_) {

    for (const auto& homologous : contig_variant.second->getVector()) {

      for (const auto &variant : homologous->getMap()) {

        variant_vector.push_back(variant.second);

      }

    }

  }

}

// Always use deep copy when modifying this object.
std::shared_ptr<kgl::GenomeVariant> kgl::GenomeVariant::deepCopy() const {

  std::shared_ptr<GenomeVariant> copy(std::make_shared<GenomeVariant>(genomeId(), ploidy()));
  copy->attributes(attributes());

  for (auto contig : getMap()) {

    std::shared_ptr<ContigVariant> copy_contig(std::make_shared<ContigVariant>(contig.second->contigId(), contig.second->ploidy()));
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
                                           PhaseId_t phase,
                                           ContigOffset_t start,
                                           ContigOffset_t end,
                                           OffsetVariantMap& variant_map) const {

  std::shared_ptr<ContigVariant> contig_variant_ptr;
  if (not getContigVariant(contig_id, contig_variant_ptr)) {

    ExecEnv::log().error("Contig Id: {} not found in Genome Variant: {}", contig_id, genomeId());
    return false;

  }

  return contig_variant_ptr->getSortedVariants(phase, start, end, variant_map);

}

/*
bool kgl::GenomeVariant::getCodingSortedVariants(std::shared_ptr<const CodingSequence> coding_sequence_ptr,
                                                 OffsetVariantMap& variant_map,
                                                 bool& frame_shift) const {

  OffsetVariantMap all_variant_map;
  if (not getSortedVariants(coding_sequence_ptr->contig()->contigId(),
                            coding_sequence_ptr->start(),
                            coding_sequence_ptr->end(),
                            all_variant_map)) {

    return false;

  }

  SignedOffset_t shift_count = 0;
  ContigOffset_t coding_sequence_offset;
  ContigSize_t coding_sequence_length;
  variant_map.clear();
  for (auto variant : all_variant_map) {

    if (SequenceOffset::refOffsetWithinCodingSequence(coding_sequence_ptr,
                                                      variant.second->offset(),
                                                      coding_sequence_offset,
                                                      coding_sequence_length)) {

      variant_map.insert(std::pair<ContigOffset_t, std::shared_ptr<const Variant>>(variant.first, variant.second));
      if (variant.second->isDelete()) {

        shift_count -= variant.second->size();

      } else if (variant.second->isInsert()) {

        shift_count += variant.second->size();

      }

    }


  }

  frame_shift = (shift_count % Codon::CODON_SIZE) != 0;

  return true;

}
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::PhasedPopulation::getGenomeVariant(const GenomeId_t& genome_id,
                                              std::shared_ptr<const GenomeVariant>& genome_variant) const {

  auto result = population_variant_map_.find(genome_id);

  if (result != population_variant_map_.end()) {

    genome_variant = result->second;
    return true;

  } else {

    genome_variant = nullptr;
    return false;

  }

}


bool kgl::PhasedPopulation::addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant) {

  auto result = population_variant_map_.insert(std::pair<GenomeId_t, std::shared_ptr<const GenomeVariant>>(genome_variant->genomeId(), genome_variant));

  return result.second;

}


size_t kgl::PhasedPopulation::variantCount() const {

  size_t variant_count = 0;
  for (auto genome : getMap()) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}


