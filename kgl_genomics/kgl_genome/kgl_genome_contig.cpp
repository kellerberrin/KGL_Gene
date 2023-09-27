//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kgl_genome_contig.h"

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


// Given a gene id and an mRNA id (sequence id) return the coding base sequence.
bool kgl::ContigReference::getCodingSequence(const FeatureIdent_t& gene_id,
                                             const FeatureIdent_t& sequence_id,
                                             std::shared_ptr<const TranscriptionSequence>& coding_sequence_ptr) const {

  std::vector<std::shared_ptr<const Feature>> feature_ptr_vec;
  std::shared_ptr<const GeneFeature> gene_ptr;
  if (findFeatureId(gene_id, feature_ptr_vec)) {

    for (const auto& feature_ptr : feature_ptr_vec) {

      if (feature_ptr->id() == gene_id) {

        if (feature_ptr->isGene()) {

          gene_ptr = std::dynamic_pointer_cast<const GeneFeature>(feature_ptr);
          break;

        } else {

          ExecEnv::log().warn("Feature: {} is not a gene.", gene_id);
          return false;

        }

      }

    }

  }

  if (not gene_ptr) {

    ExecEnv::log().warn("Gene not found for feature id: {}.", gene_id);
    return false;

  }

  auto sequence_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);
  if (sequence_array_ptr->empty()) {

    ExecEnv::log().warn("No valid coding sequences found for Gene: {}.", gene_id);
    return false;

  }

  for (const auto& sequence: sequence_array_ptr->getMap()) {

    if (sequence.second->getParent()->id() == sequence_id) {

      coding_sequence_ptr = sequence.second;
      return true;

    }

  }

  ExecEnv::log().warn("No valid coding sequences found for sequence id: {}.", sequence_id);
  return false;

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
