//
// Created by kellerberrin on 22/11/17.
//


#include "kgl_variant_factory_compound.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound SNP variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Compound SNPs are created to hold 2 or 3 coding SNPs that modify an Amino Acid simultaneously.
// 1. The SNP can only be in the same coding sequence (same mRNA parent).
// 2. They must be within an offset of 2 of each other.
// 3. They must lie within the same codon boundary (a subset of point 2).



// The logic of this function is very convoluted so read the following carefully before modifying this function.

// This function aggregates variants that occur consecutively in coding sequences (only).
// Aggregate/Compound variants are not defined for non-coding sequence variants.
// Aggregate variants exist because of the aggregate interaction with the underlying codon structure of the protein.

// 1. It is important to remember that variants in coding sequences (only) are unique to each coding sequence
// 2. There can be multiple coding sequences defined over a contig (chromosome) interval (within a gene).
// 3. Therefore there can be multiple SNP variants defined for a particular coding sequence contig offset.
// 4. These can only be distinguished by examining the coding sequence membership.
// 5. Thus coding sequences can create multiple aggregated variants for the same contig (chromosome) interval.
// 6. For SNPs there can be multiple SNPs (up to 5) per offset in a coding sequence. This is different from insert/delete.

bool kgl::CompoundSNPFactory::aggregateVariants(const std::shared_ptr<const GenomeVariant>& variant_ptr,
                                                std::vector<std::shared_ptr<const CompoundVariantMap>>& aggregated_variants_vec) const {

  // The working structure that maintains a vector of compound variants.
  std::vector<std::shared_ptr<CompoundVariantMap>> compound_variant_vec;

  for (const auto& contig_variants : variant_ptr->getMap()) {  // For all contigs,

    // flush the working structure for every contig.
    compound_variant_vec.clear();

    for (const auto& variant : contig_variants.second->getMap()) {  // For all variants.

      // only coding sequence and single SNPs (no recursive compound)
      if (variant.second->type() == VariantSequenceType::CDS_CODING and variant.second->isSingle() and variant.second->isSNP()) {

        // Get a (virtual) subordinate SNP pointer.
        std::shared_ptr<const SingleVariant> subSNP_ptr = std::static_pointer_cast<const SingleVariant>(variant.second);

          // If empty then just add the variant to a variant map and update the working structure.
          if (compound_variant_vec.empty()) {

            CompoundVariantMap compound_variant;
            std::pair<ContigSize_t, std::shared_ptr<const SingleVariant>> insert_pair(subSNP_ptr->offset(), subSNP_ptr);
            compound_variant.insert(insert_pair);
            compound_variant_vec.push_back(std::make_shared<CompoundVariantMap>(compound_variant));

          } else { // if not empty

            // scroll through all the compound variants
            bool found_sequence = false;
            std::vector<std::shared_ptr<CompoundVariantMap>> temp_compound_variant_vec;
            for (const auto &compound_variant_ptr : compound_variant_vec) {

              // Does the SNP belong to the coding sequence.
              if (compound_variant_ptr->rbegin()->second->codingSequenceId() == variant.second->codingSequenceId()) {

                // Get the codon offset of the SNP in the compoundmap and the codon offset of the SNP.

                ContigOffset_t variant_codon_offset;
                ContigOffset_t compound_codon_offset;
                ContigSize_t base_in_codon;

                if (not subSNP_ptr->codonOffset(variant_codon_offset, base_in_codon)) {

                  ExecEnv::log().error("aggregateVariants(), unexpected - could not find codon offset for variant: {}",
                                       subSNP_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
                  return false;

                }
                if (not compound_variant_ptr->rbegin()->second->codonOffset(compound_codon_offset, base_in_codon)) {

                  ExecEnv::log().error("aggregateVariants(), unexpected - could not find codon offset for variant: {}",
                                       compound_variant_ptr->rbegin()->second->output(' ', VariantOutputIndex::START_0_BASED, true));
                  return false;

                }

                if (compound_codon_offset == variant_codon_offset) { // In the same codon

                  // set the found flag to true.
                  found_sequence = true;
                  // create the insert pair.
                  std::pair<ContigSize_t, std::shared_ptr<const SingleVariant>> insert_pair(subSNP_ptr->offset(),
                                                                                             subSNP_ptr);

                  // Is the variant at the same contig address.
                  if (compound_variant_ptr->rbegin()->second->offset() == subSNP_ptr->offset()) {

                    // used if an identical offset is found. take a copy of the current map.
                    std::shared_ptr<CompoundVariantMap> temp_map_ptr(std::make_shared<CompoundVariantMap>(*compound_variant_ptr));
                    // pop the final SNP from the copy.
                    temp_map_ptr->erase(std::prev(temp_map_ptr->end()));
                    // and insert this SNP in it's place.
                    auto result = temp_map_ptr->insert(insert_pair);
                    if (not result.second) {

                      ExecEnv::log().error("aggregateVariants(), unexpected - could not add variant: {}, previous variant: {}",
                                           subSNP_ptr->output(' ', VariantOutputIndex::START_0_BASED, true),
                                           temp_map_ptr->rbegin()->second->output(' ', VariantOutputIndex::START_0_BASED, true));

                    } else { // insertion successful

                    //  ExecEnv::log().info("aggregateVariants(), inserted duplicate offset variant: {}",
                    //                       subSNP_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
                      temp_compound_variant_vec.push_back(temp_map_ptr);

                    }

                  }
                  else {

                    auto result = compound_variant_ptr->insert(insert_pair);
                    if (not result.second) {

                      ExecEnv::log().error("aggregateVariants(), unexpected - could not add variant: {}, previous variant: {}",
                                           subSNP_ptr->output(' ', VariantOutputIndex::START_0_BASED, true),
                                           compound_variant_ptr->rbegin()->second->output(' ', VariantOutputIndex::START_0_BASED, true));

                    } // could not insert

                  } // different contig offset

                } // same codon offset.

              } // if same CDS sequence

            }  // for all compound variants.

            // flush non contiguous compound variants from the working structure.
            for (const auto &compound_variant_ptr : compound_variant_vec) {

              // if the variant is not contigous with ANY compound variant then remove the compound variant.
              if ((compound_variant_ptr->rbegin()->second->contigOffset() + 2) < variant.second->contigOffset()) {

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
              std::pair<ContigSize_t, std::shared_ptr<const SingleVariant>> insert_pair(subSNP_ptr->offset(), subSNP_ptr);
              compound_variant.insert(insert_pair);
              compound_variant_vec.push_back(std::make_shared<CompoundVariantMap>(compound_variant));

            } // if not found coding sequence

          } // if not empty

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



std::shared_ptr<kgl::Variant> kgl::CompoundSNPFactory::createCompoundVariant(const CompoundVariantMap& variant_map) const {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  Phred_t quality = calculateQuality(variant_map);

  // create the variant
  std::shared_ptr<Variant> compound_insert(std::make_shared<CompoundSNP>(variant_map.begin()->second->variantSource(),
                                                                         variant_map.begin()->second->contig(),
                                                                         variant_map.begin()->second->contigOffset(),
                                                                         quality,
                                                                         variant_map));


  return compound_insert;

}




