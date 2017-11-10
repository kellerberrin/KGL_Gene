//
// Created by kellerberrin on 31/10/17.
//

#include <memory>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::ContigVariant>
kgl::ContigVariant::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::ContigVariant> filtered_contig_ptr(std::make_shared<kgl::ContigVariant>(*this));
  // Complements the bool returned by filterVariant(filter) because the delete pattern expects bool true for deletion.
  auto predicate = [&](const OffsetVariantMap::const_iterator& it) { return not it->second->filterVariant(filter); };
  predicateIterableDelete(filtered_contig_ptr->offset_variant_map_,  predicate);

  return filtered_contig_ptr;

}

// This function will insert multiple variants for each CDS sequence within each gene.
void kgl::ContigVariant::addVariant(std::shared_ptr<const Variant>& variant_ptr) {

  // Annotate the variant with genome information.
  GeneVector gene_vector;
  ContigOffset_t variant_offset = variant_ptr->contigOffset();
  if (variant_ptr->contig()->findGenes(variant_offset, gene_vector)) {

    for (const auto& gene_ptr : gene_vector) {

      std::shared_ptr<const CodingSequenceArray> sequence_array = kgl::GeneFeature::getCodingSequences(gene_ptr);
      if (sequence_array->size() == 0) {

        std::const_pointer_cast<Variant>(variant_ptr)->defineIntron(gene_ptr); // intron
        offset_variant_map_.insert(std::make_pair(variant_ptr->contigOffset(), variant_ptr));

      } else {

        for (const auto& sequence : sequence_array->getMap()) {

          if (sequence.second->isWithinCoding(variant_offset)) {

            std::const_pointer_cast<Variant>(variant_ptr)->defineCoding(sequence.second); // coding
            offset_variant_map_.insert(std::make_pair(variant_ptr->contigOffset(), variant_ptr));

          } else {  // an intron for this sequence

            std::const_pointer_cast<Variant>(variant_ptr)->defineIntron(gene_ptr); // intron
            offset_variant_map_.insert(std::make_pair(variant_ptr->contigOffset(), variant_ptr));

          } // if valid sequence for offset

        } // for all sequences within a gene

      } // if gene has a valid sequence.

    } // for all genes.

  } else {

    std::const_pointer_cast<Variant>(variant_ptr)->defineNonCoding(); // non coding
    offset_variant_map_.insert(std::make_pair(variant_ptr->contigOffset(), variant_ptr));

  }


}

std::string kgl::ContigVariant::output(char delimter) const {

  std::string variant_text;
  for (const auto& variant : offset_variant_map_) {

    variant_text += variant.second->output(delimter);

  }

  return variant_text;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GenomeVariant::addContigVariant(std::shared_ptr<kgl::ContigVariant>& contig_variant) {

  auto result = genome_variant_map_.insert(std::make_pair(contig_variant->contigId(), contig_variant));

  return result.second;

}

bool kgl::GenomeVariant::getContigVariant(const ContigId_t& contig_id,
                                          std::shared_ptr<ContigVariant>& contig_variant) {
  bool result;

  auto contig = genome_variant_map_.find(contig_id);

  if (contig != genome_variant_map_.end()) {

    contig_variant = contig->second;
    result = true;

  } else {

    contig_variant = nullptr;
    result = false;

  }

  return result;

}


bool kgl::GenomeVariant::addVariant(std::shared_ptr<const Variant> variant) {

  std::shared_ptr<ContigVariant> contig_variant;
  if (not getContigVariant(variant->contigId(), contig_variant)) {

    ExecEnv::log().error("Contig: {} not found, variant: {}", variant->contigId(), variant->output(' '));
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
kgl::GenomeVariant::emptyGenomeVariant(const VariantType_t& variant_type,
                                       const GenomeId_t& genome_id,
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




std::shared_ptr<kgl::GenomeVariant>
kgl::GenomeVariant::disaggregateCompoundVariants(const std::shared_ptr<const GenomeDatabase>& genome_db) const {

  std::shared_ptr<kgl::GenomeVariant>
  disaggreagated = emptyGenomeVariant("disaggregated variants", genomeId(), genome_db);

  for (auto contig_variant : genome_variant_map_) {

    for (auto variant : contig_variant.second->getMap()) {

      if (variant.second->isCompound()) {

        std::shared_ptr<const CompoundVariant>
        compound_variant = std::static_pointer_cast<const CompoundVariant>(variant.second);

        for (auto single_variant : compound_variant->getMap()) {

          if (not disaggreagated->addVariant(single_variant.second)) {

            ExecEnv::log().error("Cannot add disaggregated variant: {} - same contig offset as existing variant",
                                 single_variant.second->output(' '));

          }

        }

      }

    }

  }

  return disaggreagated;

}


size_t kgl::GenomeVariant::size() const {

  size_t total_variants = 0;
  for (auto contig_variant : genome_variant_map_) {

    total_variants += contig_variant.second->size();

  }

  return total_variants;

}

std::string kgl::GenomeVariant::output(char field_delimter) const {

  std::string variant_text;
  for (const auto& contig_variant : genome_variant_map_) {

    variant_text += contig_variant.second->output(field_delimter);

  }

  return variant_text;

}

std::ostream& operator<<(std::ostream &os, const kgl::GenomeVariant& genome_variant) {

  os << genome_variant.output(' ');
  os.flush();

  return os;

}
