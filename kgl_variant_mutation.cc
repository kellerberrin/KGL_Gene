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


  std::shared_ptr<DNA5SequenceCoding> DNA_reference;
  std::vector<std::shared_ptr<DNA5SequenceCoding>> DNA_mutant_vector;
  if (not mutantCodingDNA(contig_id, gene_id, sequence_id, genome_db, DNA_reference, DNA_mutant_vector)) {

    ExecEnv::log().warn("mutantProtein(), Problem generating stranded mutant DNA");
    return false;

  }

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // The reference Amino sequence.
  reference_sequence = contig_ptr->getAminoSequence(DNA_reference);

  // Mutant Amino vector
  for (const auto& mutant_sequence : DNA_mutant_vector) {

    mutant_sequence_vector.push_back(contig_ptr->getAminoSequence(mutant_sequence));

  }

  return true;

}


bool kgl::GenomeVariant::mutantCodingDNA( const ContigId_t& contig_id,
                                          const FeatureIdent_t& gene_id,
                                          const FeatureIdent_t& sequence_id,
                                          const std::shared_ptr<const GenomeDatabase>& genome_db,
                                          std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                                          std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector) const {
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
  std::shared_ptr<DNA5SequenceCoding> stranded_sequence_ptr;
  if (not contig_ptr->getDNA5SequenceCoding(coding_sequence_ptr, stranded_sequence_ptr))  {

    ExecEnv::log().warn("No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  // Return the reference amino sequence.
  reference_sequence = stranded_sequence_ptr;

  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap coding_variant_map;
  getCodingSortedVariants(contig_id, coding_sequence_ptr->start(), coding_sequence_ptr->end(), coding_variant_map);

  // There may be more than one different variant specified per offset.
  // If this is the case, then we create alternative mutation paths.
  // In other words, there can be more than one mutant protein alternative.
  // The number of mutation is exponential. Alternatives = 2 ^ (#equal offset variants) - assuming 2 alternatives per offset.
  // A message is issued if there are more than 32 alternative mutation paths (5 equal offset variants).
  // For performance reasons, a hard limit of 128 alternatives is imposed (can be varied).
  std::vector<OffsetVariantMap> variant_map_vector;
  size_t alternative_count = 0;
  std::shared_ptr<const OffsetVariantMap> variant_map_ptr(std::make_shared<OffsetVariantMap>(coding_variant_map));
  getMutationAlternatives(variant_map_ptr, variant_map_vector, alternative_count, MUTATION_SOFT_LIMIT_, MUTATION_HARD_LIMIT_);

  for (auto variant_map : variant_map_vector) {

    // Mutate.
    std::shared_ptr<DNA5SequenceCoding> mutant_coding_dna;
    if (not GenomeVariant::mutateDNA(variant_map, contig_ptr, coding_sequence_ptr,  mutant_coding_dna)) {

      ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                          contig_id, gene_id, sequence_id);
      return false;

    }

    mutant_sequence_vector.push_back(mutant_coding_dna);

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


// Returns a maximum of MUTATION_HARD_LIMIT_ alternative mutations.
bool kgl::GenomeVariant::mutantRegion( const ContigId_t& contig_id,
                                       const ContigOffset_t & region_offset,
                                       const ContigSize_t region_size,
                                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                                       std::shared_ptr<DNA5SequenceLinear>& reference_sequence,
                                       std::vector<std::shared_ptr<DNA5SequenceLinear>>& mutant_sequence_vector) const {

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the reference DNA sequence
  reference_sequence = contig_ptr->sequence().unstrandedRegion(region_offset, region_size);


  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap region_variant_map;
  getSortedVariants(contig_id, region_offset, region_offset + region_size, region_variant_map);

  // There may be more than one different variant specified per offset.
  // If this is the case, then we create alternative mutation paths.
  // In other words, there can be more than one mutant protein alternative.
  // The number of mutation is exponential. Alternatives = 2 ^ (#equal offset variants) - assuming 2 alternatives per offset.
  // A message is issued if there are more than 32 alternative mutation paths (5 equal offset variants).
  // For performance reasons, a hard limit of 128 alternatives is imposed (can be varied).
  std::vector<OffsetVariantMap> variant_map_vector;
  size_t alternative_count = 0;
  std::shared_ptr<const OffsetVariantMap> variant_map_ptr(std::make_shared<OffsetVariantMap>(region_variant_map));
  getMutationAlternatives(variant_map_ptr, variant_map_vector, alternative_count, MUTATION_SOFT_LIMIT_, MUTATION_HARD_LIMIT_);

  IndelAccountingMap indel_accounting_map;

  for (auto variant_map : variant_map_vector) {

    // Make a copy of the linear dna.
    std::shared_ptr<DNA5SequenceLinear> copy_dna_sequence_ptr(std::make_shared<DNA5SequenceLinear>(*reference_sequence));

    // And mutate it.
    if (not GenomeVariant::mutateDNA(variant_map, region_offset, copy_dna_sequence_ptr, indel_accounting_map)) {

      ExecEnv::log().warn("Problem mutating region DNA sequence for contig: {}, offset: {}, size: {}",
                          contig_id, region_offset, region_size);
      return false;

    }

    mutant_sequence_vector.push_back(copy_dna_sequence_ptr);

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

bool kgl::GenomeVariant::mutateDNA(const OffsetVariantMap& region_variant_map,
                                   ContigOffset_t contig_offset,
                                   std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr,
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


// Mutate the DNA sequence using SNP variants
bool kgl::GenomeVariant::mutateSNPs(const OffsetVariantMap& snp_variant_map,
                                    ContigOffset_t contig_offset,
                                    std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr) {

  // Mutate the base sequence.
  for (const auto& variant : snp_variant_map) {

    if (variant.second->isSingle()) {

      return mutateSingleSNP(variant.second, contig_offset, dna_sequence_ptr);

    } else { // compound SNP.

      std::shared_ptr<const CompoundSNP> cmp_snp_ptr = std::dynamic_pointer_cast<const CompoundSNP>(variant.second);

      if (not cmp_snp_ptr) {

        ExecEnv::log().error("mutateSNPs(), should be CompoundSNP; unexpected variant: {}",
                             variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      } else {

        for (auto variant : cmp_snp_ptr->getMap()) {

          if (not mutateSingleSNP(variant.second, contig_offset, dna_sequence_ptr)) {

            ExecEnv::log().error("mutateSNPs(), problem mutating sequence with compound snp: {}",
                                 cmp_snp_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
            return false;
          }

        }

      }

    }

  }

  return true;

}

// Mutate the DNA sequence using SNP variants
bool kgl::GenomeVariant::mutateSingleSNP(std::shared_ptr<const Variant> variant_ptr,
                                         ContigOffset_t contig_offset,
                                         std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr) {

  std::shared_ptr<const SNPVariant> snp_ptr = std::dynamic_pointer_cast<const SNPVariant>(variant_ptr);

  if (not snp_ptr) {

    ExecEnv::log().error("mutateSingleSNP, should be SNPVariant; unexpected variant: {}",
                         variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }

  // Adjust the offset
  SignedOffset_t adjusted_offset = snp_ptr->offset() - contig_offset;

  // Check the offset
  if (adjusted_offset < 0 or adjusted_offset >= dna_sequence_ptr->length()) {

    ExecEnv::log().error("mutateSingleSNP(), calculated sequence offset: {} is out of range for sequence size: {}",
                         adjusted_offset, dna_sequence_ptr->length(),
                         snp_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  ContigOffset_t sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);

  // Check the reference.
  if (snp_ptr->reference() != dna_sequence_ptr->at(sequence_offset)) {

    ExecEnv::log().warn("mutateSingleSNP(), reference base: {} does not match sequence base: {} at sequence offset: {}",
                        DNA5::convertToChar(snp_ptr->reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        sequence_offset);

  }

  // Mutate the sequence
  dna_sequence_ptr->modifyBase(sequence_offset, snp_ptr->mutant());

  return true;

}

// Mutate the DNA sequence using Delete variants
bool kgl::GenomeVariant::mutateDeletes(const OffsetVariantMap& snp_variant_map,
                                       ContigOffset_t contig_offset,
                                       std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr,
                                       IndelAccountingMap& indel_accounting_map) {

  return true;

}


// Mutate the DNA sequence using Insert variants
bool kgl::GenomeVariant::mutateInserts(const OffsetVariantMap& snp_variant_map,
                                       ContigOffset_t contig_offset,
                                       std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr,
                                       IndelAccountingMap& indel_accounting_map) {

  return true;

}
