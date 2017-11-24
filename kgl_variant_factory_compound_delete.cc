//
// Created by kellerberrin on 22/11/17.
//


#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound delete variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::VariantDeleteFactory::selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const {

  if (variant_ptr->isSNP()) {

    const std::shared_ptr<const SNPVariantDNA5> SNP_ptr = std::static_pointer_cast<const SNPVariantDNA5>(variant_ptr);
    return ExtendDNA5::isDeletion(SNP_ptr->mutant());

  } else {

    return false;

  }

}



std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantDeleteFactory::compoundDelete(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                          const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<kgl::GenomeVariant>
  compound_delete_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_variants->genomeId(), genome_db_ptr);

  std::vector<std::shared_ptr<const CompoundVariantMap>> aggregated_variants_vec;

  // Get the aggregated variants.
  aggregateVariants(genome_variants ,aggregated_variants_vec);

  // Generate the actual compound deletes.
  generateCompoundDeletes(aggregated_variants_vec, compound_delete_variants);

  return compound_delete_variants;

}


void kgl::VariantDeleteFactory::generateCompoundDeletes( const std::vector<std::shared_ptr<const CompoundVariantMap>>& contiguous_delete_vec,
                                                         std::shared_ptr<kgl::GenomeVariant>& genome_variant_ptr) {

  for (const auto& variant_map : contiguous_delete_vec) {

    std::shared_ptr<const Variant> delete_variant = createCompoundDelete(*variant_map);
    if (delete_variant != nullptr) {

      if (not genome_variant_ptr->addVariant(delete_variant)) {

        ExecEnv::log().error("Unable to add compound delete variant: {} - probable offset duplicate",
                             delete_variant->output(' ', VariantOutputIndex::START_0_BASED));

      }

    }

  } // for all contiguous

}


std::shared_ptr<const kgl::Variant>
kgl::VariantDeleteFactory::createCompoundDelete(const CompoundVariantMap& variant_map) {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  if (variant_map.begin()->second->geneMembership().empty()) {

    ExecEnv::log().error("Variant has no associated gene");
    return nullptr;

  }

  std::shared_ptr<const GeneFeature> gene_ptr = variant_map.begin()->second->geneMembership().front();

  StrandSense gene_strand = gene_ptr->sequence().strand();
  std::shared_ptr<const ContigFeatures> contig_ptr;
  ContigOffset_t variant_offset;

  switch(gene_strand) {

    case StrandSense::UNKNOWN:
    case StrandSense::FORWARD:
      contig_ptr = variant_map.begin()->second->contig();
      variant_offset = variant_map.begin()->second->contigOffset();
      break;

    case StrandSense::REVERSE:
      contig_ptr = variant_map.rbegin()->second->contig();
      variant_offset = variant_map.rbegin()->second->contigOffset();
      break;

    default: // cannot happen, but the compiler complains.
      contig_ptr = variant_map.begin()->second->contig();
      variant_offset = variant_map.begin()->second->contigOffset();
      break;

  }

  // source is the same as the snp variants
  const std::string& variant_source = variant_map.begin()->second->variantSource();
  // create the variant
  std::shared_ptr<Variant> compound_delete(std::make_shared<CompoundDelete>(variant_source,
                                                                            contig_ptr,
                                                                            variant_offset,
                                                                            variant_map));
  // define its coding sequence.
  compound_delete->defineCoding(variant_map.begin()->second->codingSequences().getFirst());

  return compound_delete;

}
