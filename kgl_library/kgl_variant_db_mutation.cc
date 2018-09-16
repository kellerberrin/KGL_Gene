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
                                         const std::shared_ptr<const GenomeDatabase>& genome_db,
                                         OffsetVariantMap& variant_map,
                                         std::shared_ptr<AminoSequence>& reference_sequence,
                                         std::shared_ptr<AminoSequence>& mutant_sequence) const {


  std::shared_ptr<DNA5SequenceCoding> DNA_reference;
  std::shared_ptr<DNA5SequenceCoding> DNA_mutant;
  if (not mutantCodingDNA(contig_id, phase, gene_id, sequence_id, genome_db, variant_map, DNA_reference, DNA_mutant)) {

    ExecEnv::log().warn("mutantProtein(), Problem generating stranded mutant DNA");
    return false;

  }

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Check the reference sequence for good measure.
  if (not DNA_reference->verifySequence()) {

    ExecEnv::log().error("mutantProteins(), corrupted reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // The reference Amino sequence.
  reference_sequence = contig_ptr->getAminoSequence(DNA_reference);

  // Check that the DNA corruption is acquired by mutating the reference sequence
  // NOT from a problem with the codon arithmetic (tested separately).
  if (not DNA_mutant->verifySequence()) {

    ExecEnv::log().error("mutantProteins(), corrupted mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);
    for (auto variant : variant_map) {

      ExecEnv::log().info("mutantProteins(), variant map: {}", variant.second->output(' ', VariantOutputIndex::START_1_BASED, true));

    }


  }

  mutant_sequence = contig_ptr->getAminoSequence(DNA_mutant);

  return true;

}


bool kgl::GenomeVariant::mutantCodingDNA( const ContigId_t& contig_id,
                                          PhaseId_t phase,
                                          const FeatureIdent_t& gene_id,
                                          const FeatureIdent_t& sequence_id,
                                          const std::shared_ptr<const GenomeDatabase>& genome_db,
                                          OffsetVariantMap& variant_map,
                                          std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                                          std::shared_ptr<DNA5SequenceCoding>& mutant_sequence) const {
  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantCodingDNA(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the coding sequence.
  std::shared_ptr<const CodingSequence> coding_sequence_ptr;
  if (not contig_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("mutantCodingDNA(), Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }


  // Get the reference DNA sequence
  std::shared_ptr<DNA5SequenceCoding> stranded_sequence_ptr;
  if (not contig_ptr->getDNA5SequenceCoding(coding_sequence_ptr, stranded_sequence_ptr))  {

    ExecEnv::log().warn("No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  // Return the reference sequence.
  reference_sequence = stranded_sequence_ptr;

  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap coding_variant_map;
  getSortedVariants(contig_id, phase, coding_sequence_ptr->start(), coding_sequence_ptr->end(), coding_variant_map);
  variant_map = coding_variant_map;

  if (not VariantMutation().mutateDNA(coding_variant_map, contig_ptr, coding_sequence_ptr, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                        contig_id, gene_id, sequence_id);
    return false;

  }

  return true;

}

// Returns a maximum of MUTATION_HARD_LIMIT_ alternative mutations.
bool kgl::GenomeVariant::mutantRegion( const ContigId_t& contig_id,
                                       PhaseId_t phase,
                                       ContigOffset_t region_offset,
                                       ContigSize_t region_size,
                                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                                       OffsetVariantMap& variant_map,
                                       std::shared_ptr<DNA5SequenceLinear>& reference_sequence,
                                       std::shared_ptr<DNA5SequenceLinear>& mutant_sequence) const {

  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantRegion(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Check offset and size.
  if ((region_offset + region_size) > contig_ptr->sequence().length() or region_size > contig_ptr->sequence().length()) {

    ExecEnv::log().warn("mutantRegion(), contig offset: {} and region size: {} too large for contig: {} length: {}",
                         region_offset, region_size, contig_ptr->contigId(), contig_ptr->sequence().length());
    return false;

  }

  // Get the reference DNA sequence
  reference_sequence = contig_ptr->sequence().unstrandedRegion(region_offset, region_size);


  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap region_variant_map;
  getSortedVariants(contig_id, phase, region_offset, region_offset + region_size, region_variant_map);
  variant_map = region_variant_map;

  // Make a copy of the linear dna.
  mutant_sequence = std::make_shared<DNA5SequenceLinear>(*reference_sequence);

  // And mutate it.
  if (not VariantMutation().mutateDNA(region_variant_map, region_offset, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating region DNA sequence for contig: {}, offset: {}, size: {}",
                        contig_id, region_offset, region_size);
    return false;

  }

  return true;

}


bool kgl::GenomeVariant::mutantContig( const ContigId_t& contig_id,
                                       PhaseId_t phase,
                                       const std::shared_ptr<const GenomeDatabase>& genome_db,
                                       std::shared_ptr<const DNA5SequenceContig>& reference_contig,
                                       std::shared_ptr<DNA5SequenceContig>& mutant_contig_ptr) const {


  // Get the contig.
  std::shared_ptr<const ContigFeatures> contig_ptr;
  if (not genome_db->getContigSequence(contig_id, contig_ptr)) {

    ExecEnv::log().warn("mutantProtein(), Could not find contig: {} in genome database", contig_id);
    return false;

  }

  // Get the contig DNA sequence
  reference_contig = contig_ptr->sequence_ptr();


  // Generate the mutant sequences.
  // Extract the variants for processing.
  OffsetVariantMap contig_variant_map;
  getSortedVariants(contig_id, phase,  0, contig_ptr->contigSize(), contig_variant_map);

  // There are no alternative paths used in mutating a contig.
  // Make a copy of the contig dna.
  mutant_contig_ptr = std::make_shared<DNA5SequenceContig>(*reference_contig);

    // And mutate it.
  if (not VariantMutation().mutateDNA(contig_variant_map, 0, mutant_contig_ptr)) {

      ExecEnv::log().warn("Problem mutating genome: {},  contig: {}", contig_id);
      return false;

  }

  return true;


}

