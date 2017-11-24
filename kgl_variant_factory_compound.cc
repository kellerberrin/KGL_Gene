//
// Created by kellerberrin on 22/11/17.
//



#include "kgl_variant_factory_compound.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<kgl::GenomeVariant>
kgl::VariantCompoundFactory::disaggregateCompoundVariants(const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                          const std::shared_ptr<const GenomeDatabase>& genome_db) const {

  std::shared_ptr<kgl::GenomeVariant> disaggreagated = genome_variant->emptyGenomeVariant(genome_variant->genomeId(),
                                                                                          genome_db);

  for (auto contig_variant : genome_variant->contigMap()) {

    for (auto variant : contig_variant.second->getMap()) {

      if (variant.second->isCompound()) {

        std::shared_ptr<const CompoundVariant>
        compound_variant = std::static_pointer_cast<const CompoundVariant>(variant.second);

        for (auto single_variant : compound_variant->getMap()) {

          if (not disaggreagated->addVariant(single_variant.second)) {

            ExecEnv::log().error("Cannot add disaggregated variant: {} - same contig offset as existing variant",
                                 single_variant.second->output(' ', VariantOutputIndex::START_0_BASED));

          }

        }

      }

    }

  }

  return disaggreagated;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Aggregate insert/delete variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The logic of this function is very convoluted so read the following carefully before modifying this function.

// This function aggregates variants that occur consecutively in coding sequences (only).
// Aggregate/Compound variants are not defined for non-coding sequence variants.
// Aggregate variants exist because of the aggregate interaction with the underlying codon structure of the protein.

// 1. It is important to remember that variants in coding sequences (only) are unique to each coding sequence
// 2. There can be multiple coding sequences defined over a contig (chromosome) interval (within a gene).
// 3. Therefore there can be multiple SNP variants defined for a particular coding sequence contig offset.
// 4. These can only be distinguished by examining the coding sequence membership.
// 5. Thus coding sequences can create multiple aggregated variants for the same contig (chromosome) interval.
bool kgl::VariantInsertDeleteFactory::aggregateVariants(const std::shared_ptr<const GenomeVariant>& variant_ptr,
                                                        std::vector<std::shared_ptr<const CompoundVariantMap>>& aggregated_variants_vec) const {

  // The working structure that maintains a vector of compound variants.
  std::vector<std::shared_ptr<CompoundVariantMap>> compound_variant_vec;

  for (const auto& contig_variants : variant_ptr->contigMap()) {  // For all contigs,

    // flush the working structure for every contig.
    compound_variant_vec.clear();

    for (const auto& variant : contig_variants.second->getMap()) {  // For all variants.

      // only coding sequence and non-compound (no recursive compound)
      if (variant.second->type() == VariantSequenceType::CDS_CODING and not variant.second->isCompound()) {

        if (selectVariant(variant.second)) {  // The selection function passed in as a template.

          // If empty then just add the variant to a variant map and update the working structure.
          if (compound_variant_vec.empty()) {

            CompoundVariantMap compound_variant;
            compound_variant.insert(variant);
            compound_variant_vec.push_back(std::make_shared<CompoundVariantMap>(compound_variant));

          } else { // if not empty

            // scroll through all the compound variants
            bool found_sequence = false;
            for (const auto& compound_variant_ptr : compound_variant_vec) {

              // Does the variant belong to the coding sequence.
              if (compound_variant_ptr->rbegin()->second->codingSequenceId() == variant.second->codingSequenceId()) {

                // set the found flag to true.
                found_sequence = true;
                // If the variant is the at same offset, then this is an error.
                if (compound_variant_ptr->rbegin()->second->contigOffset() == variant.second->contigOffset()) {

                  // Error if in the same coding sequence.
                  ExecEnv::log().error("Variant: {} with the same offset in the same coding sequence as variant: {}",
                                       compound_variant_ptr->rbegin()->second->output(' ', VariantOutputIndex::START_0_BASED),
                                       variant.second->output(' ', VariantOutputIndex::START_0_BASED));

                } else if ((compound_variant_ptr->rbegin()->second->contigOffset() + 1) == variant.second->contigOffset()) {

                  // Is contiguous in the same coding sequence so append to the compound map.
                  auto result = compound_variant_ptr->insert(variant);
                  if (not result.second) {

                    ExecEnv::log().error("aggregateCompoundVariants(), unexpected - could not insert to contiguous map");
                    return false;

                  }

                }

                break; // found the coding sequence so no need to look further.

              } // If in the same sequence.

            }  // for all compound variants.

            // flush non contiguous compound variants from the working structure.
            std::vector<std::shared_ptr<CompoundVariantMap>> temp_compound_variant_vec;
            for (const auto& compound_variant_ptr : compound_variant_vec) {

              // if the variant is not contigous with ANY compound variant then remove the compound variant.
              if ((compound_variant_ptr->rbegin()->second->contigOffset() + 1) < variant.second->contigOffset()) {

                // candidate compound variants must be at least size 2.
                if (compound_variant_ptr->size() >= 2) {

                  aggregated_variants_vec.push_back(compound_variant_ptr);

                }

              } else { // still active so copy to the temp variable.

                temp_compound_variant_vec.push_back(compound_variant_ptr);

              }

            }  // for all compound variants
            // copy the temp structure to the working structure.
            compound_variant_vec = temp_compound_variant_vec;
            // Finally, if we did not find a coding sequence for the variant, then add it to the working structure.
            if (not found_sequence) {

              CompoundVariantMap compound_variant;
              compound_variant.insert(variant);
              compound_variant_vec.push_back(std::make_shared<CompoundVariantMap>(compound_variant));

            } // if not found coding sequence

          } // if not empty

        } // if selected variant type.

      } // if in a coding sequence and non-compound (no recursive compound)

    } // for all variants

    // we are moving to another contig so flush the working structure.

    for (const auto& compound_variant_ptr : compound_variant_vec) {

      // candidate compound variants must be at least size 2.
      if (compound_variant_ptr->size() >= 2) {

        aggregated_variants_vec.push_back(compound_variant_ptr);

      } // if size > 2

    }  // for all compound variants

  } // for all contigs.

  return true;

}


