//
// Created by kellerberrin on 22/12/17.
//

#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::mutantProteins( const ContigId_t& contig_id,
                                         const FeatureIdent_t& gene_id,
                                         const FeatureIdent_t& sequence_id,
                                         const std::shared_ptr<const GenomeDatabase>& genome_db,
                                         std::shared_ptr<AminoSequence>& reference_sequence,
                                         std::vector<std::shared_ptr<AminoSequence>>& mutant_sequence_vector) const {
  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }


  // Get the reference DNA sequence
  std::shared_ptr<DNA5SequenceCoding> dna_sequence_ptr;
  if (not contig_ptr->getDNA5SequenceCoding(coding_sequence_ptr, dna_sequence_ptr))  {

    ExecEnv::log().warn("No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  // Return the reference amino sequence.
  reference_sequence = contig_ptr->getAminoSequence(dna_sequence_ptr);

  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap coding_variant_map;
  getCodingSortedVariants(contig_id, coding_sequence_ptr->start(), coding_sequence_ptr->end(), coding_variant_map);

  // There maybe more than one variant specified per offset.
  // If this is the case then we create alternative mutation paths.
  // In other words, there can be more than one mutant protein alternative.
  // This function is exponential. Alternatives = 2 ^ (#equal offset variants) - assuming 2 alternatives.
  // A message is issued if there are more than 32 alternatives (5 equal offset variants).
  // For performance reasons, a hard limit of 128 alternatives is imposed (can be varied).
  std::vector<OffsetVariantMap> variant_map_vector;
  size_t alternative_count = 0;
  std::shared_ptr<const OffsetVariantMap> variant_map_ptr(std::make_shared<OffsetVariantMap>(coding_variant_map));
  getMutationAlternatives(variant_map_ptr, variant_map_vector, alternative_count, MUTATION_SOFT_LIMIT_, MUTATION_HARD_LIMIT_);

  for (auto variant_map : variant_map_vector) {

    // Make a copy of the protein dna.
    std::shared_ptr<DNA5SequenceCoding> copy_dna_sequence_ptr(std::make_shared<DNA5SequenceCoding>(*dna_sequence_ptr));

    // And mutate it.
    if (not GenomeVariant::mutateDNA(variant_map, sequence_id, copy_dna_sequence_ptr)) {

      ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                          contig_id, gene_id, sequence_id);
      return false;

    }

    mutant_sequence_vector.push_back(contig_ptr->getAminoSequence(copy_dna_sequence_ptr));

  }

  return true;

}

// Recursively generates mutation alternatives.
// Warning - not a long function, but the coding logic is convoluted.
// Study carefully before modification.
void kgl::GenomeVariant::getMutationAlternatives(std::shared_ptr<const OffsetVariantMap> variant_map_ptr,
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


// Perform the actual mutation.
bool kgl::GenomeVariant::mutateDNA(const OffsetVariantMap& variant_map,
                                   const FeatureIdent_t& sequence_id,
                                   std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {

  // Split the variant map into SNP, Delete and Insert Variants.

  OffsetVariantMap snp_variant_map;
  OffsetVariantMap delete_variant_map;
  OffsetVariantMap insert_variant_map;
  SplitVariantMap(variant_map, snp_variant_map, delete_variant_map, insert_variant_map);

  DeleteAccountingMap delete_accounting_map;

  if (not mutateSNPs(snp_variant_map, sequence_id, dna_sequence_ptr)) {

    ExecEnv::log().error("Problem with SNP mutations");
    return false;

  }

  if (not mutateDeletes(delete_variant_map, sequence_id, dna_sequence_ptr, delete_accounting_map)) {

    ExecEnv::log().error("Problem with Delete mutations");
    return false;

  }

  if (not mutateInserts(insert_variant_map, sequence_id, delete_accounting_map, dna_sequence_ptr)) {

    ExecEnv::log().error("Problem with Insert mutations");
    return false;

  }

  return true;

}



// Split the variant map into SNP, Delete and Insert Variants.
void kgl::GenomeVariant::SplitVariantMap(const OffsetVariantMap& variant_map,
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


// Mutate the DNA sequence using SNP variants
bool kgl::GenomeVariant::mutateSNPs(const OffsetVariantMap& snp_variant_map,
                                    const FeatureIdent_t& sequence_id,
                                    std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {

  // Mutate the base sequence.
  for (const auto& variant : snp_variant_map) {

    SignedOffset_t sequence_size_adjust = 0;  // How the variant modifies sequence size.
    if (not variant.second->mutateCodingSequence(sequence_id, 0, dna_sequence_ptr->length(), sequence_size_adjust, dna_sequence_ptr)) {

      ExecEnv::log().error("mutateDNA(), problem with variant: {}",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    if (sequence_size_adjust != 0) {

      ExecEnv::log().error("mutateDNA(), Non zero sequence resize: {} with SNP variant: {}",
                           sequence_size_adjust, variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

  }

  return true;

}


// Mutate the DNA sequence using delete variants
bool kgl::GenomeVariant::mutateDeletes(const OffsetVariantMap& delete_variant_map,
                                       const FeatureIdent_t& sequence_id,
                                       std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr,
                                       DeleteAccountingMap& delete_accounting_map) {

  // Process the deletions in reverse order so that we don't have to adjust the variant offsets.

  for (auto rit = delete_variant_map.rbegin(); rit != delete_variant_map.rend(); ++rit) {



  }

  return true;

}


// Mutate the DNA sequence using insert variants
bool kgl::GenomeVariant::mutateInserts(const OffsetVariantMap& insert_variant_map,
                                       const FeatureIdent_t& sequence_id,
                                       const DeleteAccountingMap& delete_accounting_map,
                                       std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {


  return true;

}

