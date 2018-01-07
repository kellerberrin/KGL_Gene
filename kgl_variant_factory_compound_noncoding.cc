//
// Created by kellerberrin on 6/01/18.
//

#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound non coding insert variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// The logic of this function is very convoluted so read the following carefully before modifying this function.



bool kgl::NonCodingInsertDeleteFactory::aggregateVariants(const std::shared_ptr<const GenomeVariant>& genome_variant_ptr,
                                                          std::vector<std::shared_ptr<const CompoundVariantMap>>& aggregated_variants_vec) const {

  // The working structure that maintains a vector of compound variants.
  std::vector<std::shared_ptr<CompoundVariantMap>> compound_variant_vec;

  for (const auto& contig_variants : genome_variant_ptr->getMap()) {  // For all contigs,

    // flush the working structure for every contig.
    compound_variant_vec.clear();

    for (const auto& variant : contig_variants.second->getMap()) {  // For all variants.

      // polymorphic selection.
      if (selectVariant(variant.second) and variant.second->isSingle() and variant.second->type() != VariantSequenceType::CDS_CODING) {  // The variant selection

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
                ExecEnv::log().warn("Identical Variant: {} with the same offset as variant: {}",
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


void kgl::NonCodingInsertDeleteFactory::nonCodingIntron(std::shared_ptr<Variant> variant_ptr) const {

  GeneVector gene_vector;
  ContigOffset_t variant_offset = variant_ptr->contigOffset();
  if (variant_ptr->contig()->findGenes(variant_offset, gene_vector)) {

    if (not gene_vector.empty()) {

      variant_ptr->defineIntron(gene_vector.front()); // intron

      if (gene_vector.size() > 1) {

        ExecEnv::log().error("Intron Non Coding Compound Insert is a Member of Multiple genes : {}",
                             variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));

      }

    } else {

      // tag as non-coding.
      variant_ptr->defineNonCoding(); // non coding

    }

  } else {

    // tag as non-coding.
    variant_ptr->defineNonCoding(); // non coding


  }


}


std::shared_ptr<kgl::Variant> kgl::CompoundNonCodingInsertFactory::createCompoundVariant(const CompoundVariantMap& variant_map) const {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  Phred_t quality = calculateQuality(variant_map);

  // create the variant
  std::shared_ptr<Variant> compound_insert(std::make_shared<CompoundInsert>(variant_map.begin()->second->variantSource(),
                                                                         variant_map.begin()->second->contig(),
                                                                         variant_map.begin()->second->contigOffset(),
                                                                         quality,
                                                                         variant_map));

  nonCodingIntron(compound_insert);

  return compound_insert;

}


std::shared_ptr<kgl::Variant> kgl::CompoundNonCodingDeleteFactory::createCompoundVariant(const CompoundVariantMap& variant_map) const {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  Phred_t quality = calculateQuality(variant_map);

  // create the variant
  std::shared_ptr<Variant> compound_delete(std::make_shared<CompoundDelete>(variant_map.begin()->second->variantSource(),
                                                                            variant_map.begin()->second->contig(),
                                                                            variant_map.begin()->second->contigOffset(),
                                                                            quality,
                                                                            variant_map));

  nonCodingIntron(compound_delete);

  return compound_delete;

}


