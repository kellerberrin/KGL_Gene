//
// Created by kellerberrin on 3/01/18.
//

#include <memory>
#include "kgl_variant_db.h"
#include "kgl_variant_mutation.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::mutantProteins( const ContigId_t& contig_id,
                                         const FeatureIdent_t& gene_id,
                                         const FeatureIdent_t& sequence_id,
                                         const std::shared_ptr<const GenomeDatabase>& genome_db,
                                         bool& frame_shift_flag,
                                         std::shared_ptr<AminoSequence>& reference_sequence,
                                         std::vector<std::shared_ptr<AminoSequence>>& mutant_sequence_vector) const {


  std::shared_ptr<DNA5SequenceCoding> DNA_reference;
  std::vector<std::shared_ptr<DNA5SequenceCoding>> DNA_mutant_vector;
  if (not mutantCodingDNA(contig_id, gene_id, sequence_id, genome_db, frame_shift_flag, DNA_reference, DNA_mutant_vector)) {

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
                                          bool& frame_shift_flag,
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
  getCodingSortedVariants(contig_id,
                          gene_id,
                          sequence_id,
                          coding_sequence_ptr->start(),
                          coding_sequence_ptr->end(),
                          coding_variant_map,
                          frame_shift_flag);

  // There may be more than one different variant specified per offset.
  // If this is the case, then we create alternative mutation paths.
  // In other words, there can be more than one mutant protein alternative.
  // The number of mutation is exponential. Alternatives = 2 ^ (#equal offset variants) - assuming 2 alternatives per offset.
  // A message is issued if there are more than 32 alternative mutation paths (5 equal offset variants).
  // For performance reasons, a hard limit of 128 alternatives is imposed (can be varied).
  std::vector<OffsetVariantMap> variant_map_vector;
  size_t alternative_count = 0;
  std::shared_ptr<const OffsetVariantMap> variant_map_ptr(std::make_shared<OffsetVariantMap>(coding_variant_map));
  VariantMutation::getMutationAlternatives(variant_map_ptr, variant_map_vector, alternative_count);

  for (auto variant_map : variant_map_vector) {

    // Mutate if non-empty variant map.
    if (not variant_map.empty()) {

      std::shared_ptr<DNA5SequenceCoding> mutant_coding_dna;
      VariantMutation variant_mutation;
      if (not variant_mutation.mutateDNA(variant_map, contig_ptr, coding_sequence_ptr, mutant_coding_dna)) {

        ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                            contig_id, gene_id, sequence_id);
        return false;

      }

      mutant_sequence_vector.push_back(mutant_coding_dna);

    }

  }

  return true;

}

// Returns a maximum of MUTATION_HARD_LIMIT_ alternative mutations.
bool kgl::GenomeVariant::mutantRegion( const ContigId_t& contig_id,
                                       ContigOffset_t region_offset,
                                       ContigSize_t region_size,
                                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                                       std::shared_ptr<DNA5SequenceLinear>& reference_sequence,
                                       std::vector<std::shared_ptr<DNA5SequenceLinear>>& mutant_sequence_vector) const {

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Check offset and size.
  if ((region_offset + region_size) >= contig_ptr->sequence().length() or region_size >= contig_ptr->sequence().length()) {

    ExecEnv::log().warn("mutantRegion(), contig offset: {} and region size: {} too large for contig: {} length: {}",
                         region_offset, region_size, contig_ptr->contigId(), contig_ptr->sequence().length());
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
  VariantMutation::getMutationAlternatives(variant_map_ptr, variant_map_vector, alternative_count);

  for (auto variant_map : variant_map_vector) {

    // Non-empty variant map.
    if (not variant_map.empty()) {

      // Make a copy of the linear dna.
      std::shared_ptr<DNA5SequenceLinear> copy_dna_sequence_ptr(std::make_shared<DNA5SequenceLinear>(*reference_sequence));

      // And mutate it.
      VariantMutation variant_mutation;
      if (not variant_mutation.mutateDNA(variant_map, region_offset, copy_dna_sequence_ptr)) {

        ExecEnv::log().warn("Problem mutating region DNA sequence for contig: {}, offset: {}, size: {}",
                            contig_id, region_offset, region_size);
        return false;

      }

      mutant_sequence_vector.push_back(copy_dna_sequence_ptr);

    }

  }

  return true;

}

