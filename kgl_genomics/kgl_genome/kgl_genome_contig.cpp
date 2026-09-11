//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kgl_genome_contig.h"

#include <ranges>
#include <optional>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contig Reference members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::ContigReference::ContigReference(ContigId_t contig_id,
                                      const std::shared_ptr<const DNA5SequenceLinear>& sequence_ptr)
                                      : contig_id_(std::move(contig_id)), sequence_ptr_(sequence_ptr) {

  if (not sequence_ptr_) {

    ExecEnv::log().critical("ContigReference::ContigReference; contig id: {} has a null sequence pointer", contig_id_);

  }

}


bool kgl::ContigReference::addContigFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  return gene_exon_features_.checkAddFeature(feature_ptr);

}


void kgl::ContigReference::verifyFeatureHierarchy() {

  // Setup the Gene feature structure first.
  gene_exon_features_.setupVerifyHierarchy();
  // Verify the Genes.
  verifyGeneFeatures();

}


// Convenience routine for Amino sequences.
kgl::AminoSequence kgl::ContigReference::getAminoSequence(const DNA5SequenceCoding& coding_sequence) const {

  return coding_table_.getAminoSequence(coding_sequence);

}


std::vector<std::shared_ptr<const kgl::Feature>> kgl::ContigReference::findFeatureId(const FeatureIdent_t& feature_id) const {

  return gene_exon_features_.findFeatureId(feature_id);

}


// Given a gene id and an mRNA id (sequence id) return the coding base sequence.
std::optional<std::shared_ptr<const kgl::TranscriptionSequence>>
kgl::ContigReference::getTranscription(const FeatureIdent_t& gene_id, const FeatureIdent_t& transcript_id) const {

  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec = findFeatureId(gene_id);

  auto find_iter = std::ranges::find_if(feature_ptr_vec, [](const std::shared_ptr<const Feature>& feature_ptr) { return feature_ptr->isGene(); });
  if (find_iter == feature_ptr_vec.end()) {

    ExecEnv::log().warn("Feature id: {} is not a Gene", gene_id);
    return std::nullopt;

  }

  auto gene_ptr = std::dynamic_pointer_cast<const GeneFeature>(*find_iter);
  if (not gene_ptr) {

    ExecEnv::log().error("ContigReference::getTranscription; Feature id: {} is marked as a gene but is not a GeneFeature", gene_id);
    return std::nullopt;

  }

  auto transcript_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);
  for (const auto& transcript_ptr : transcript_array_ptr->getMap() | std::views::values) {

    if (transcript_ptr->getParent()->id() == transcript_id) {

      return transcript_ptr;

    }

  }

  ExecEnv::log().warn("No valid coding sequences found for Gene: {}, transcript id: {}.", gene_id, transcript_id);
  return std::nullopt;

}


bool kgl::ContigReference::equivalent(const ContigReference& lhs) const {

  return contig_id_ == lhs.contig_id_
         and *sequence_ptr_ == *lhs.sequence_ptr_
         and gene_exon_features_.equivalent(lhs.gene_exon_features_)
         and coding_table_.translationTableName() == lhs.coding_table_.translationTableName();

}

// Given a gene transcript, generate the associated (strand adjusted) coding sequence.
std::optional<kgl::DNA5SequenceCoding>
kgl::ContigReference::codingSequence(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) const {

  const auto cds_interval_set = transcript_ptr->getExonIntervals();
  auto concat_sequence_opt = sequence().concatSequences(cds_interval_set);
  if (not concat_sequence_opt) {

    ExecEnv::log().warn("Unable to concat sequence intervals for Gene: {}, Transcript: {}",
                        transcript_ptr->getGene()->id(), transcript_ptr->getParent()->id());

    return std::nullopt; // Return an empty coding sequence.

  }
  const auto& concat_sequence = concat_sequence_opt.value();

  return concat_sequence.codingSequence(transcript_ptr->strand());

}
