//
// Created by kellerberrin on 23/04/18.
//


#include <memory>
#include "kel_patterns.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;



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

