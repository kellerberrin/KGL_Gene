//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  if (feature_ptr->isTSS()) {

    aux_contig_features_.checkAddFeature(feature_ptr);

  } else {

    gene_exon_features_.checkAddFeature(feature_ptr);

  }

  return true;

}



void kgl::ContigFeatures::verifyFeatureHierarchy() {

  // Setup the Gene feature structure first.
  gene_exon_features_.setupVerifyHierarchy();
  // Verify the Genes.
  verifyCDSPhasePeptide();
  // Setup the Aux feature hierarchy using the gene exon hierarchy.
  aux_contig_features_.setupVerifyHierarchy(gene_exon_features_);

}



// Convenience routine for Amino sequences.
std::shared_ptr<kgl::AminoSequence>
kgl::ContigFeatures::getAminoSequence(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

  return coding_table_.getAminoSequence(sequence_ptr);

}


// Given a gene id and an mRNA id (sequence id) return the coding base sequence.
bool kgl::ContigFeatures::getCodingSequence(const FeatureIdent_t& gene_id,
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
bool kgl::ContigFeatures::getDNA5SequenceCoding(const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                                                std::shared_ptr<DNA5SequenceCoding>& sequence_ptr) const {

  if (coding_sequence_ptr) {

    sequence_ptr = sequence_ptr_->DNA5SequenceContig::codingSequence(coding_sequence_ptr);
    return true;

  }

  ExecEnv::log().error("getDNA5SequenceCoding(), coding_sequence_ptr is null");
  return false;

}

