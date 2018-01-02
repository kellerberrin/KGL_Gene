//
// Created by kellerberrin on 22/12/17.
//

#include <memory>
#include "kgl_variant_single.h"
#include "kgl_variant_compound.h"
#include "kgl_sequence_virtual_compare.h"
#include "kgl_variant_mutation.h"

namespace kgl = kellerberrin::genome;


// Mutate coding sequences
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& variant_map,
                                     std::shared_ptr<const ContigFeatures> contig_ptr,
                                     std::shared_ptr<const CodingSequence> coding_sequence_ptr,
                                     std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {

  // Split the variant map into SNP, Delete and Insert Variants.

  OffsetVariantMap snp_variant_map;
  OffsetVariantMap delete_variant_map;
  OffsetVariantMap insert_variant_map;
  SplitVariantMap(variant_map, snp_variant_map, delete_variant_map, insert_variant_map);

  IndelAccountingMap indel_accounting_map;

  // mutate UNSTRANDED DNA and then convert to STRANDED DNA
  std::shared_ptr<DNA5SequenceLinear>
  unstranded_ptr = contig_ptr->sequence().unstrandedRegion(coding_sequence_ptr->start(),
                                                           (coding_sequence_ptr->end() - coding_sequence_ptr->start()));

  if (not mutateSNPs(snp_variant_map, coding_sequence_ptr->start(), unstranded_ptr)) {

    ExecEnv::log().error("Problem with SNP mutations");
    return false;

  }

  if (not mutateDeletes(delete_variant_map, coding_sequence_ptr->start(), unstranded_ptr, indel_accounting_map)) {

    ExecEnv::log().error("Problem with Delete mutations");
    return false;

  }

  if (not mutateInserts(insert_variant_map, coding_sequence_ptr->start(),  unstranded_ptr, indel_accounting_map)) {

    ExecEnv::log().error("Problem with Insert mutations");
    return false;

  }

  dna_sequence_ptr = unstranded_ptr->codingSubSequence(coding_sequence_ptr, 0, 0, coding_sequence_ptr->start());

  return true;

}


// Mutate DNA region.
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& region_variant_map,
                                     ContigOffset_t contig_offset,
                                     std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                     IndelAccountingMap& indel_accounting_map) {

  // Split the variant map into SNP, Delete and Insert Variants.

  OffsetVariantMap snp_variant_map;
  OffsetVariantMap delete_variant_map;
  OffsetVariantMap insert_variant_map;
  SplitVariantMap(region_variant_map, snp_variant_map, delete_variant_map, insert_variant_map);

  indel_accounting_map.clear();

  if (not mutateSNPs(snp_variant_map, contig_offset, dna_sequence_ptr)) {

    ExecEnv::log().error("Problem with SNP mutations");
    return false;

  }

  if (not mutateDeletes(delete_variant_map, contig_offset, dna_sequence_ptr, indel_accounting_map)) {

    ExecEnv::log().error("Problem with Delete mutations");
    return false;

  }

  if (not mutateInserts(insert_variant_map, contig_offset, dna_sequence_ptr, indel_accounting_map)) {

    ExecEnv::log().error("Problem with Insert mutations");
    return false;

  }

  return true;

}


// Recursively generates mutation alternatives.
// Warning - not a long function, but the coding logic is convoluted.
// Study carefully before modification.
void kgl::VariantMutation::getMutationAlternatives(std::shared_ptr<const OffsetVariantMap> variant_map_ptr,
                                                   std::vector<OffsetVariantMap>& variant_map_vector,
                                                   size_t& alternative_count,
                                                   size_t soft_limit,
                                                   size_t hard_limit) {

  alternative_count += 1;

  if (alternative_count == soft_limit) {

    ExecEnv::log().info("Soft limit of {} coding sequence mutations reached", soft_limit);

  }

  if (alternative_count > hard_limit) {

    ExecEnv::log().warn("Hard limit of {} coding sequence mutations reached, no additional mutations generated", hard_limit);
    return;

  }

  OffsetVariantMap alternative_map;
  for (auto it = variant_map_ptr->begin(); it != variant_map_ptr->end(); ++it) {

    if (not alternative_map.empty()) {

      if (alternative_map.rbegin()->second->offsetOverlap(*it->second)) {

        std::shared_ptr<OffsetVariantMap> copy_map_ptr(std::make_shared<OffsetVariantMap>(alternative_map));
        copy_map_ptr->erase(std::prev(copy_map_ptr->end()));  // Pop the last element

        for (auto copy_it = it; copy_it != variant_map_ptr->end(); ++copy_it) { // Copy the rest.

          copy_map_ptr->insert(*copy_it);

        }
        // Recursive call.
        getMutationAlternatives(copy_map_ptr, variant_map_vector, alternative_count, soft_limit, hard_limit);

      } else { // no matching offset variant.

        alternative_map.insert(*it);

      }

    } else { // empty

      alternative_map.insert(*it);

    }

  }

  variant_map_vector.push_back(alternative_map);

}


// Split the variant map into SNP, Delete and Insert Variants.
void kgl::VariantMutation::SplitVariantMap(const OffsetVariantMap& variant_map,
                                         OffsetVariantMap& snp_variant_map,
                                         OffsetVariantMap& delete_variant_map,
                                         OffsetVariantMap& insert_variant_map) {

  snp_variant_map.clear();
  delete_variant_map.clear();
  insert_variant_map.clear();

  for (auto variant : variant_map) {

    if (variant.second->isSNP()) {

      snp_variant_map.insert(variant);

    } else if (variant.second->isDelete()) {

      delete_variant_map.insert(variant);

    } else if (variant.second->isInsert()) {

      insert_variant_map.insert(variant);

    } else { // Unknown variant type.

      ExecEnv::log().error("SplitVariantMap(), Unknown variant type :{}",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return;
    }

  }

}

