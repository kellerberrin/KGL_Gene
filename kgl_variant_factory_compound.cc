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
kgl::CompoundFactory::disaggregate(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                          const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) const {

  std::shared_ptr<kgl::GenomeVariant> disaggreagated = genome_variants->emptyGenomeVariant(genome_variants->genomeId(),
                                                                                           genome_db_ptr);

  for (auto contig_variant : genome_variants->getMap()) {

    for (auto variant : contig_variant.second->getMap()) {

      if (variant.second->isCompound()) {

        std::shared_ptr<const CompoundVariant>
        compound_variant = std::static_pointer_cast<const CompoundVariant>(variant.second);

        for (auto single_variant : compound_variant->getMap()) {

          disaggreagated->addVariant(single_variant.second);

        }

      }

    }

  }

  return disaggreagated;

}




std::shared_ptr<kgl::GenomeVariant>
kgl::CompoundFactory::create(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                    const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) const {

  std::shared_ptr<kgl::GenomeVariant>
  compound_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_variants->genomeId(), genome_db_ptr);

  std::vector<std::shared_ptr<const CompoundVariantMap>> aggregated_variants_vec;

// Get the aggregated variants.
  aggregateVariants(genome_variants, aggregated_variants_vec);

// Add to the genome variant pointer.
  for (const auto &variant_map : aggregated_variants_vec) {

    std::shared_ptr<Variant> compound_variant_ptr = createCompoundVariant(*variant_map);

    if (compound_variant_ptr) {

      if (not compound_variants->addVariant(compound_variant_ptr)) {

        ExecEnv::log().error("Unable to add compound variant: {} - probable offset duplicate",
                             compound_variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));

      }

    }

  } // for all variant maps

  return compound_variants;

}


// Assumes the component subordinate quality probabilities are independent.
// Needs to be re-coded.
kgl::Phred_t kgl::CompoundFactory::calculateQuality(const CompoundVariantMap& variant_map) const
{

  double probability = 0.0;
  double adjust_log;
  for (auto variant : variant_map) {

    adjust_log = variant.second->quality() / -10.0;
    probability += ::pow10(adjust_log);

  }

  adjust_log = ::log10(probability) * -10.0;

  return adjust_log > 0 ? adjust_log : 0.0;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Generate the compound variant maps for insert/delete
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

bool kgl::InsertDeleteFactory::aggregateVariants(const std::shared_ptr<const GenomeVariant>& genome_variant_ptr,
                                                 std::vector<std::shared_ptr<const CompoundVariantMap>>& aggregated_variants_vec) const {

  // The working structure that maintains a vector of compound variants.
  std::vector<std::shared_ptr<CompoundVariantMap>> compound_variant_vec;

  for (const auto& contig_variants : genome_variant_ptr->getMap()) {  // For all contigs,

    // flush the working structure for every contig.
    compound_variant_vec.clear();

    for (const auto& variant : contig_variants.second->getMap()) {  // For all variants.

      // only coding sequence and non-compound (no recursive compound)
      if (variant.second->type() == VariantSequenceType::CDS_CODING and variant.second->isSingle()) {

        if (selectVariant(variant.second)) {  // The virtual selection function.

          std::shared_ptr<const SingleVariant> subvariant_ptr = std::dynamic_pointer_cast<const SingleVariant>(variant.second);

          if (not subvariant_ptr) {

            ExecEnv::log().error("Not a (Subordinate) Variant: {}",
                                 variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
            return false;

          }

          // If empty then just add the variant to a variant map and update the working structure.
          if (compound_variant_vec.empty()) {

            std::pair<ContigOffset_t, std::shared_ptr<const SingleVariant>> insert_pair(subvariant_ptr->offset(), subvariant_ptr);
            CompoundVariantMap compound_variant;
            compound_variant.insert(insert_pair);
            compound_variant_vec.push_back(std::make_shared<CompoundVariantMap>(compound_variant));

          } else { // if not empty

            // scroll through all the compound variants
            bool found_sequence = false;
            std::vector<std::shared_ptr<CompoundVariantMap>> add_compound_variant_vec;
            for (const auto& compound_variant_ptr : compound_variant_vec) {

              // Does the variant belong to the coding sequence.
              if (compound_variant_ptr->rbegin()->second->codingSequenceId() == variant.second->codingSequenceId()) {

                // If the variant is the at same offset
                if (compound_variant_ptr->rbegin()->second->contigOffset() == variant.second->contigOffset()) {

                  // If it is a different variant then copy the compound variant map, pop the last entry, push the variant and add the map.
                  if (not compound_variant_ptr->rbegin()->second->equivalent(*variant.second)) {

                    CompoundVariantMap compound_variant = *compound_variant_ptr; // copy
                    compound_variant.erase(std::prev(compound_variant.end())); // pop last
                    std::pair<ContigSize_t, std::shared_ptr<const SingleVariant>> insert_pair(subvariant_ptr->offset(), subvariant_ptr);
                    compound_variant.insert(insert_pair); // push this
                    add_compound_variant_vec.push_back(std::make_shared<CompoundVariantMap>(compound_variant)); // add to working structure.

                  } else {

                    // Error if an identical variant in the same coding sequence.
                    ExecEnv::log().warn("Identical Variant: {} with the same offset in the same coding sequence as variant: {}",
                                        compound_variant_ptr->rbegin()->second->output(' ', VariantOutputIndex::START_0_BASED, true),
                                        variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));

                  }

                } else if ((compound_variant_ptr->rbegin()->second->contigOffset() + 1) == variant.second->contigOffset()) {

                  // set the found flag to true.
                  found_sequence = true;
                  // Is contiguous in the same coding sequence so append to the compound map.
                  std::pair<ContigSize_t, std::shared_ptr<const SingleVariant>> insert_pair(subvariant_ptr->offset(), subvariant_ptr);
                  auto result = compound_variant_ptr->insert(insert_pair);
                  if (not result.second) {

                    ExecEnv::log().error("aggregateCompoundVariants(), unexpected - could not insert to contiguous map");
                    return false;

                  }

                }

              } // If in the same sequence.

            }  // for all compound variants.

            // add in any additional compound variants.
            for (const auto& compound_variant : add_compound_variant_vec) {

              compound_variant_vec.push_back(compound_variant);

            }

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
              std::pair<ContigSize_t, std::shared_ptr<const SingleVariant>> insert_pair(subvariant_ptr->offset(), subvariant_ptr);
              compound_variant.insert(insert_pair);
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


