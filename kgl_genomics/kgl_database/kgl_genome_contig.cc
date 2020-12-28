//
// Created by kellerberrin on 12/11/17.
//


#include "kel_exec_env.h"
#include "kel_patterns.h"
#include "kgl_genome_contig.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigReference::addGeneExonFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  return gene_exon_features_.checkAddFeature(feature_ptr);

}


bool kgl::ContigReference::addAuxFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  aux_contig_features_.checkAddFeature(feature_ptr);

  return true;

}




void kgl::ContigReference::verifyFeatureHierarchy() {

  // Setup the Gene feature structure first.
  gene_exon_features_.setupVerifyHierarchy();
  // Verify the Genes.
  verifyCDSPhasePeptide();

}


void kgl::ContigReference::verifyAuxillaryHierarchy() {

  // Setup the Aux feature hierarchy using the gene exon hierarchy.
  aux_contig_features_.setupVerifyHierarchy(gene_exon_features_);

}


// Convenience routine for Amino sequences.
kgl::AminoSequence kgl::ContigReference::getAminoSequence(const DNA5SequenceCoding& coding_sequence) const {

  return coding_table_.getAminoSequence(coding_sequence);

}


// Given a gene id and an mRNA id (sequence id) return the coding base sequence.
bool kgl::ContigReference::getCodingSequence(const FeatureIdent_t& gene_id,
                                             const FeatureIdent_t& sequence_id,
                                             std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const {

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

  std::shared_ptr<const CodingSequenceArray> sequence_array_ptr = GeneFeature::getCodingSequences(gene_ptr);

  if (sequence_array_ptr->empty()) {

    ExecEnv::log().warn("No valid coding sequences found for Gene: {}.", gene_id);
    return false;

  }

  for (const auto& sequence: sequence_array_ptr->getMap()) {

    if (sequence.second->getCDSParent()->id() == sequence_id) {

      coding_sequence_ptr = sequence.second;
      return true;

    }

  }

  ExecEnv::log().warn("No valid coding sequences found for sequence id: {}.", sequence_id);
  return false;

}


// Given a CDS coding sequence, return the corresponding DNA base sequence (strand adjusted).
bool kgl::ContigReference::getDNA5SequenceCoding(const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                                                 DNA5SequenceCoding& coding_sequence) const {

  if (coding_sequence_ptr) {

    coding_sequence = sequence_ptr_->DNA5SequenceContig::codingSequence(coding_sequence_ptr);
    return true;

  }

  ExecEnv::log().error("getDNA5SequenceCoding(), coding_sequence_ptr is null");
  return false;

}

