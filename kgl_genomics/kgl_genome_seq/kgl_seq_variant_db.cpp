//
// Created by kellerberrin on 3/01/18.
//

#include <memory>
#include "kgl_seq_variant.h"
#include "kgl_seq_coding.h"
#include "kgl_seq_variant_db.h"


namespace kgl = kellerberrin::genome;


bool kgl::GenomeMutation::mutantProteins( const ContigId_t& contig_id,
                                         const FeatureIdent_t& gene_id,
                                         const FeatureIdent_t& sequence_id,
                                         const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                         const OffsetVariantMap& variant_map,
                                         AminoSequence& reference_sequence,
                                         AminoSequence& mutant_sequence) {


  DNA5SequenceCoding DNA_reference;
  DNA5SequenceCoding DNA_mutant;
  if (not mutantCodingDNA(contig_id, gene_id, sequence_id, genome_ref_ptr, variant_map, DNA_reference, DNA_mutant)) {

    ExecEnv::log().warn("Problem generating stranded mutant DNA");
    return false;

  }

  // Get the contig_ref_ptr.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("Could not find contig_ref_ptr: {} in genome database", contig_id);
    return false;

  }
  auto contig_ref_ptr = contig_ref_opt.value();

  // Check the DNA5 reference sequence.
  if (not DNA_reference.verifySequence()) {

    ExecEnv::log().error("Corrupted CodingDNA5 reference sequence: {}, gene: {}, contig_ref_ptr: {}", sequence_id, gene_id, contig_id);

  }

  // The reference Amino sequence.
  reference_sequence = contig_ref_ptr->getAminoSequence(DNA_reference);

  // Check the reference amino sequence.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("Corrupted amino reference sequence: {}, gene: {}, contig_ref_ptr: {}", sequence_id, gene_id, contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not DNA_mutant.verifySequence()) {

    ExecEnv::log().error("Corrupted mutant sequence: {}, gene: {}, contig_ref_ptr: {}", sequence_id, gene_id, contig_id);
    // List out the variants.
    for (auto const& variant : variant_map) {

      ExecEnv::log().info("Variant map: {}", variant.second->HGVS_Phase());

    }

  }

  mutant_sequence = contig_ref_ptr->getAminoSequence(DNA_mutant);

  // Check the mutant amino sequence.
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("Corrupted amino mutant sequence: {}, gene: {}, contig_ref_ptr: {}", sequence_id, gene_id, contig_id);

  }

  return true;

}


bool kgl::GenomeMutation::mutantCodingDNA( const ContigId_t& contig_id,
                                          const FeatureIdent_t& gene_id,
                                          const FeatureIdent_t& transcript_id,
                                          const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                          const OffsetVariantMap& variant_map,
                                          DNA5SequenceCoding& reference_sequence,
                                          DNA5SequenceCoding& mutant_sequence) {
  // Get the contig_ref_ptr.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("Could not find contig_ref_ptr: {} in genome database", contig_id);
    return false;

  }
  auto& contig_ref_ptr = contig_ref_opt.value();

  // Get the coding sequence.
  auto transcript_opt = contig_ref_ptr->getCodingSequence(gene_id, transcript_id);
  if (not transcript_opt) {

    ExecEnv::log().warn("Could not find a coding sequence for gene: {}, sequence: {}", gene_id, transcript_id);
    return false;

  }
  auto& transcript_ptr = transcript_opt.value();

  auto coding_dna_opt = contig_ref_ptr->codingSequence(transcript_ptr);
  if (not coding_dna_opt)  {

    ExecEnv::log().warn("No valid DNA sequence for contig_ref_ptr: {}, gene: {}, sequence id: {}", contig_id, gene_id, transcript_id);
    return false;

  }
  reference_sequence = std::move(coding_dna_opt.value());

  if (not VariantMutation().mutateDNA(variant_map, contig_ref_ptr, transcript_ptr, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating DNA sequence for contig_ref_ptr: {}, gene: {}, sequence id: {}", contig_id, gene_id, transcript_id);
    return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("Corrupted reference sequence: {}, gene: {}, contig_ref_ptr: {}", transcript_id, gene_id, contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("Corrupted mutant sequence: {}, gene: {}, contig_ref_ptr: {}", transcript_id, gene_id, contig_id);
    // List out the variants.
    for (auto const& [offset, variant_ptr] : variant_map) {

      ExecEnv::log().info("Variant map: {}", variant_ptr->HGVS_Phase());

    }

  }

  return true;

}

bool kgl::GenomeMutation::mutantRegion( const ContigId_t& contig_id,
                                       ContigOffset_t region_offset,
                                       ContigSize_t region_size,
                                       const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                       const OffsetVariantMap& variant_map,
                                       DNA5SequenceLinear& reference_sequence,
                                       DNA5SequenceLinear& mutant_sequence) {

  // Get the contig_ref_ptr.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("Could not find contig_ref_ptr: {} in genome database", contig_id);
    return false;

  }

  auto contig_ref_ptr = contig_ref_opt.value();

  // Check offset and size.
  if ((region_offset + region_size) > contig_ref_ptr->sequence().length() or region_size > contig_ref_ptr->sequence().length()) {

    ExecEnv::log().warn("Contig offset: {} and region size: {} too large for contig_ref_ptr: {} length: {}",
                        region_offset, region_size, contig_ref_ptr->contigId(), contig_ref_ptr->sequence().length());
    return false;

  }

  // Get the reference DNA sequence
  OpenRightUnsigned reference_interval(region_offset, region_offset+region_size);
  auto reference_sequence_opt = contig_ref_ptr->sequence().subSequence(reference_interval);
  if (not reference_sequence_opt) {

    ExecEnv::log().warn("Problem extracting DNA sequence interval:{} for contig_ref_ptr: {}, interval: {}",
                        reference_interval.toString(), contig_id, contig_ref_ptr->sequence_ptr()->interval().toString());
    return false;

  }
  reference_sequence = std::move(reference_sequence_opt.value());
  // And mutate the sequence.
  if (not VariantMutation().mutateDNA(variant_map, contig_ref_ptr, region_offset, region_size, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating region DNA sequence for contig_ref_ptr: {}, interval: {}",
                        contig_id, reference_interval.toString());
    return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("Corrupted reference contig_ref_ptr: {}, offset: {}, size: {}", contig_id, region_offset, region_size);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("Corrupted mutant contig_ref_ptr: {}, offset: {}, size: {}", contig_id, region_offset, region_size);
    // List out the variants.
    for (auto const& [offset, variant_ptr] : variant_map) {

      ExecEnv::log().info("Variant map: {}", variant_ptr->HGVS_Phase());

    }

  }

  return true;

}

