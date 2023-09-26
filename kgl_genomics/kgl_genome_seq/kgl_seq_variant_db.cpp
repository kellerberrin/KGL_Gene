//
// Created by kellerberrin on 3/01/18.
//

#include <memory>
#include "kgl_genome_seq/kgl_seq_variant.h"
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

    ExecEnv::log().warn("GenomeMutation::mutantProteins; Problem generating stranded mutant DNA");
    return false;

  }

  // Get the contig.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("GenomeMutation::mutantProteins; Could not find contig: {} in genome database", contig_id);
    return false;

  }
  auto contig_ref_ptr = contig_ref_opt.value();

  // Check the DNA5 reference sequence.
  if (not DNA_reference.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantProteins; corrupted CodingDNA5 reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // The reference Amino sequence.
  reference_sequence = contig_ref_ptr->getAminoSequence(DNA_reference);

  // Check the reference amino sequence.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantProteins; corrupted amino reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not DNA_mutant.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantProteins; corrupted mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);
    // List out the variants.
    for (auto const& variant : variant_map) {

      ExecEnv::log().info("GenomeMutation::mutantProteins; variant map: {}", variant.second->HGVS_Phase());

    }

  }

  mutant_sequence = contig_ref_ptr->getAminoSequence(DNA_mutant);

  // Check the mutant amino sequence.
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantProteins; corrupted amino mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  return true;

}


bool kgl::GenomeMutation::mutantCodingDNA( const ContigId_t& contig_id,
                                          const FeatureIdent_t& gene_id,
                                          const FeatureIdent_t& sequence_id,
                                          const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                          const OffsetVariantMap& variant_map,
                                          DNA5SequenceCoding& reference_sequence,
                                          DNA5SequenceCoding& mutant_sequence) {
  // Get the contig.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("GenomeMutation::mutantCodingDNA; Could not find contig: {} in genome database", contig_id);
    return false;

  }

  auto contig_ref_ptr = contig_ref_opt.value();

  // Get the coding sequence.
  std::shared_ptr<const TranscriptionSequence> coding_sequence_ptr;
  if (not contig_ref_ptr->getCodingSequence(gene_id, sequence_id, coding_sequence_ptr)) {

    ExecEnv::log().warn("GenomeMutation::mutantCodingDNA; Could not find a coding sequence for gene: {}, sequence: {}", gene_id, sequence_id);
    return false;

  }


  if (not contig_ref_ptr->getDNA5SequenceCoding(coding_sequence_ptr, reference_sequence))  {

    ExecEnv::log().warn("GenomeMutation::mutantCodingDNA; No valid DNA sequence for contig: {}, gene: {}, sequence id: {}", contig_id, gene_id, sequence_id);
    return false;

  }

  if (not VariantMutation().mutateDNA(variant_map, contig_ref_ptr, coding_sequence_ptr, mutant_sequence)) {

    ExecEnv::log().warn("GenomeMutation::mutantCodingDNA; Problem mutating DNA sequence for contig: {}, gene: {}, sequence id: {}",
                        contig_id, gene_id, sequence_id);
    return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantCodingDNA; corrupted reference sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantCodingDNA; corrupted mutant sequence: {}, gene: {}, contig: {}", sequence_id, gene_id, contig_id);
    // List out the variants.
    for (auto const& [offset, variant_ptr] : variant_map) {

      ExecEnv::log().info("GenomeMutation::mutantCodingDNA; variant map: {}", variant_ptr->HGVS_Phase());

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

  // Get the contig.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("GenomeMutation::mutantRegion; Could not find contig: {} in genome database", contig_id);
    return false;

  }

  auto contig_ref_ptr = contig_ref_opt.value();

  // Check offset and size.
  if ((region_offset + region_size) > contig_ref_ptr->sequence().length() or region_size > contig_ref_ptr->sequence().length()) {

    ExecEnv::log().warn("GenomeMutation::mutantRegion; contig offset: {} and region size: {} too large for contig: {} length: {}",
                        region_offset, region_size, contig_ref_ptr->contigId(), contig_ref_ptr->sequence().length());
    return false;

  }

  // Get the reference DNA sequence
  OpenRightUnsigned reference_interval(region_offset, region_offset+region_size);
  auto reference_sequence_opt = contig_ref_ptr->sequence().subOptSequence(reference_interval);
  if (not reference_sequence_opt) {

    ExecEnv::log().warn("Problem extracting DNA sequence interval:{} for contig: {}, interval: {}",
                        reference_interval.toString(), contig_id, contig_ref_ptr->sequence_ptr()->interval().toString());
    return false;

  }
  reference_sequence = std::move(reference_sequence_opt.value());
  // And mutate the sequence.
  if (not VariantMutation().mutateDNA(variant_map, contig_ref_ptr, region_offset, region_size, mutant_sequence)) {

    ExecEnv::log().warn("Problem mutating region DNA sequence for contig: {}, interval: {}",
                        contig_id, reference_interval.toString());
    return false;

  }

  // Check the reference sequence for good measure.
  if (not reference_sequence.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantRegion; corrupted reference contig: {}, offset: {}, size: {}",
                         contig_id, region_offset, region_size);

  }

  // Check for any DNA corruption acquired by mutating the reference sequence
  if (not mutant_sequence.verifySequence()) {

    ExecEnv::log().error("GenomeMutation::mutantRegion; corrupted mutant contig: {}, offset: {}, size: {}", contig_id, region_offset, region_size);
    // List out the variants.
    for (auto const& [offset, variant_ptr] : variant_map) {

      ExecEnv::log().info("GenomeMutation::mutantRegion; variant map: {}", variant_ptr->HGVS_Phase());

    }

  }

  return true;

}

