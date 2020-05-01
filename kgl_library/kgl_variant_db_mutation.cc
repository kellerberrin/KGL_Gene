//
// Created by kellerberrin on 3/01/18.
//

#include <memory>
#include "kgl_variant_db.h"
#include "kgl_variant_mutation.h"

namespace kgl = kellerberrin::genome;


bool kgl::GenomeVariant::mutantProteins( const ContigId_t& contig_id,
                                         PhaseId_t phase,
                                         const FeatureIdent_t& gene_id,
                                         const FeatureIdent_t& sequence_id,
                                         const std::shared_ptr<const GenomeReference>& genome_db,
                                         OffsetVariantMap& variant_map,
                                         AminoSequence& reference_sequence,
                                         AminoSequence& mutant_sequence) const {


  DNA5SequenceCoding DNA_reference;
  DNA5SequenceCoding DNA_mutant;
  if (not mutantCodingDNA(contig_id, phase, gene_id, sequence_id, genome_db, variant_map, DNA_reference, DNA_mutant)) {

    ExecEnv::log().warn("mutantProtein(), Problem generating stranded mutant DNA");
    return false;

  }

  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Check the DNA5 reference sequence.
  if (not DNA_reference.verifySequence()) {

    ExecEnv::log().error("mutantProteins(), corrupted CodingDNA5 reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // The reference Amino sequence.
  reference_sequence = contig_opt.value()->getAminoSequence(DNA_reference);

  // Check the reference amino sequence.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("mutantProteins(), corrupted amino reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not DNA_mutant.verifySequence()) {

    ExecEnv::log().error("mutantProteins(), corrupted mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);
    // List out the variants.
    for (auto const& variant : variant_map) {

      ExecEnv::log().info("mutantProteins(), variant map: {}", variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));

    }

  }

  mutant_sequence = contig_opt.value()->getAminoSequence(DNA_mutant);

  // Check the mutant amino sequence.
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("mutantProteins(), corrupted amino mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  return true;

}


bool kgl::GenomeVariant::mutantCodingDNA( const ContigId_t& contig_id,
                                          PhaseId_t phase,
                                          const FeatureIdent_t& gene_id,
                                          const FeatureIdent_t& sequence_id,
                                          const std::shared_ptr<const GenomeReference>& genome_db,
                                          OffsetVariantMap& variant_map,
                                          DNA5SequenceCoding& reference_sequence,
                                          DNA5SequenceCoding& mutant_sequence) const {
  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("mutantCodingDNA(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_opt.value()->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("mutantCodingDNA(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }


  if (not contig_opt.value()->getDNA5SequenceCoding(coding_sequence_ptr, reference_sequence))  {

    ExecEnv::log().warn("No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap coding_variant_map;

  if (not getSortedVariants( contig_id,
                             phase,
                             coding_sequence_ptr->start(),
                             coding_sequence_ptr->end(),
                             coding_variant_map)) {

    ExecEnv::log().error("mutantCodingDNA(), could not get sorted variants for, gene: {}, contig: {}", gene_id, contig_id);

  }

  variant_map = coding_variant_map;

  if (not VariantMutation().mutateDNA(coding_variant_map, contig_opt.value(), coding_sequence_ptr, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                        contig_id, gene_id, sequence_id);
    return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("mutantCodingDNA(), corrupted reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("mutantCodingDNA(), corrupted mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);
    // List out the variants.
    for (auto const& variant : variant_map) {

      ExecEnv::log().info("mutantCodingDNA(), variant map: {}", variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));

    }

  }

  return true;

}

bool kgl::GenomeVariant::mutantRegion( const ContigId_t& contig_id,
                                       PhaseId_t phase,
                                       ContigOffset_t region_offset,
                                       ContigSize_t region_size,
                                       const std::shared_ptr<const GenomeReference>& genome_db,
                                       OffsetVariantMap& variant_map,
                                       DNA5SequenceLinear& reference_sequence,
                                       DNA5SequenceLinear& mutant_sequence) const {

  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("mutantRegion(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Check offset and size.
  if ((region_offset + region_size) > contig_opt.value()->sequence().length() or region_size > contig_opt.value()->sequence().length()) {

    ExecEnv::log().warn("mutantRegion(), contig offset: {} and region size: {} too large for contig: {} length: {}",
                        region_offset, region_size, contig_opt.value()->contigId(), contig_opt.value()->sequence().length());
    return false;

  }

  // Get the reference DNA sequence
  reference_sequence = contig_opt.value()->sequence().subSequence(region_offset, region_size);

  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap region_variant_map;

  if (not getSortedVariants(contig_id, phase, region_offset, region_offset + region_size, region_variant_map)) {

    ExecEnv::log().warn("Problem retrieving sorted variants to mutate region DNA sequence for contig: {}, offset: {}, size: {}",
                        contig_id, region_offset, region_size);
    return false;

  }

  variant_map = region_variant_map;

  // And mutate the sequence.
  if (not VariantMutation().mutateDNA(region_variant_map, contig_opt.value(), region_offset, region_size, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating region DNA sequence for contig: {}, offset: {}, size: {}",
                        contig_id, region_offset, region_size);
    return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("mutantRegion(), corrupted reference contig: {}, offset: {}, size: {}",
                         contig_id, region_offset, region_size);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("mutantRegion(), corrupted mutant contig: {}, offset: {}, size: {}", contig_id, region_offset, region_size);
    // List out the variants.
    for (auto const& variant : variant_map) {

      ExecEnv::log().info("mutantRegion(), variant map: {}", variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));

    }

  }

  return true;

}


bool kgl::GenomeVariant::mutantContig( const ContigId_t& contig_id,
                                       PhaseId_t phase,
                                       const std::shared_ptr<const GenomeReference>& genome_db,
                                       std::shared_ptr<const DNA5SequenceContig>& reference_contig_ptr,
                                       DNA5SequenceContig& mutant_contig) const {


  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the contig DNA sequence
  reference_contig_ptr = contig_opt.value()->sequence_ptr();


  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap contig_variant_map;

  if (not getSortedVariants(contig_id, phase, 0, contig_opt.value()->contigSize(), contig_variant_map)) {

    ExecEnv::log().warn("Problem retrieving sorted variants to mutate region DNA sequence contig: {}", contig_id);
    return false;

  }

    // And mutate it.
  if (not VariantMutation().mutateDNA(contig_variant_map, contig_opt.value(), mutant_contig)) {

      ExecEnv::log().warn("Problem mutating genome: {},  contig: {}", contig_id);
      return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_contig_ptr->verifySequence()) {

    ExecEnv::log().error("mutantContig(), corrupted reference contig: {}", contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference contig
  if (not mutant_contig.verifySequence()) {

    ExecEnv::log().error("mutantContig(), corrupted mutant contig: {}", contig_id);

  }

  return true;


}

