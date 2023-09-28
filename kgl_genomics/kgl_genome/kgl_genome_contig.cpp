//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kgl_genome_contig.h"
#include "kgl_seq_interval.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigReference::addContigFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  return gene_exon_features_.checkAddFeature(feature_ptr);

}



void kgl::ContigReference::verifyFeatureHierarchy() {

  // Setup the Gene feature structure first.
  gene_exon_features_.setupVerifyHierarchy();
  // Verify the Genes.
  verifyCDSPhasePeptide();

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
  kgl::ContigReference::getCodingSequence(const FeatureIdent_t& gene_id, const FeatureIdent_t& sequence_id) const {

  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec = findFeatureId(gene_id);
  std::shared_ptr<const GeneFeature> gene_ptr;
  if (not feature_ptr_vec.empty()) {

    for (const auto& feature_ptr : feature_ptr_vec) {

      if (feature_ptr->id() == gene_id) {

        if (feature_ptr->isGene()) {

          gene_ptr = std::dynamic_pointer_cast<const GeneFeature>(feature_ptr);
          break;

        } else {

          ExecEnv::log().warn("Feature: {} is not a gene.", gene_id);
          return std::nullopt;

        }

      }

    }

  }

  if (not gene_ptr) {

    ExecEnv::log().warn("Gene not found for feature id: {}.", gene_id);
    return std::nullopt;

  }

  auto transcript_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);
  if (transcript_array_ptr->empty()) {

    ExecEnv::log().warn("No valid coding sequences found for Gene: {}.", gene_id);
    return std::nullopt;

  }

  for (const auto& [transcript_id, transcript_ptr] : transcript_array_ptr->getMap()) {

    if (transcript_ptr->getParent()->id() == sequence_id) {

      return transcript_ptr;

    }

  }

  ExecEnv::log().warn("No valid coding sequences found for sequence id: {}.", sequence_id);
  return std::nullopt;

}


bool kgl::ContigReference::equivalent(const ContigReference& lhs) const {

  if (contig_id_ != lhs.contig_id_) {

    return false;

  }

  bool compare_sequence = *sequence_ptr_ == *lhs.sequence_ptr_;
  if (not compare_sequence) {

    return false;

  }

  if (not gene_exon_features_.equivalent(lhs.gene_exon_features_)) {

    return false;

  }

  if (coding_table_.translationTableName() != lhs.coding_table_.translationTableName()) {

    return false;

  }

  return true;

}


std::optional<kgl::DNA5SequenceCoding>
kgl::ContigReference::codingSequence( const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) {

  auto cds_interval_set = GeneIntervalStructure::transcriptIntervals(transcript_ptr);
  std::vector<OpenRightUnsigned> interval_vector(cds_interval_set.begin(), cds_interval_set.end());

  auto concat_sequence_opt = sequence().concatSequences(interval_vector);
  if (not concat_sequence_opt) {

    ExecEnv::log().warn("Unable to concat sequence intervals for Gene: {}, Transcript: {}",
                        transcript_ptr->getGene()->id(), transcript_ptr->getParent()->id());

    return std::nullopt; // Return an empty coding sequence.

  }
  auto& concat_sequence = concat_sequence_opt.value();

  return concat_sequence.codingSequence(transcript_ptr->strand());

}


kgl::ProteinSequenceValidity
kgl::ContigReference::checkValidCodingSequence(const DNA5SequenceCoding& coding_sequence) const {


  auto sequence_length = coding_sequence.length();
  if ((sequence_length % Codon::CODON_SIZE) != 0) {

    return ProteinSequenceValidity::NOT_MOD3;

  }

  auto amino_sequence = getAminoSequence(coding_sequence);
  return checkValidProteinSequence(amino_sequence);

}

