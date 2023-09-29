//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_genome_feature.h"
#include "kgl_genome_contig.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::ContigReference::verifyCDSPhasePeptide() {

  // Iterate through all the features looking for Genes.
  size_t protein_gene_count{0};
  size_t ill_formed_genes{0};
  size_t ncRNA_gene_count{0};
  size_t empty_genes{0};
  size_t valid_sequence{0};
  size_t empty_sequence{0};
  size_t not_mod3{0};
  size_t no_start_codon{0};
  size_t no_stop_codon{0};
  size_t nonsense_mutation{0};
  bool verbose{false};

  if (verbose) {

    ExecEnv::log().info("ContigReference::verifyCDSPhasePeptide; Verifying {} Genes using amino translation table: {}", contigId(), coding_table_.translationTableName());

  }

  for(const auto& [feature_id, feature_ptr  ] : gene_exon_features_.offsetFeatureMap()) {

    if (feature_ptr->isGene()) {

      auto gene_ptr = std::static_pointer_cast<const GeneFeature>(feature_ptr);
      auto sequence_array_ptr = kgl::GeneFeature::getTranscriptionSequences(gene_ptr);

      if (sequence_array_ptr->empty()) {

        ++empty_genes;
        continue;

      }

      bool valid_gene{false};
      bool ncRNA_gene{false};
      for (auto const& [parent_id, sequence_ptr]: sequence_array_ptr->getMap()) {

        if (sequence_ptr->codingType() == TranscriptionSequenceType::NCRNA) {

          ++ncRNA_gene_count;
          ncRNA_gene = true;
          continue;

        }

        switch(TranscriptionSequence::checkValidProtein(sequence_ptr)) {

          case ProteinSequenceValidity::VALID:
            valid_gene = true; // Only require 1 valid sequence for a valid protein gege.
            ++valid_sequence;
            break;

          case ProteinSequenceValidity::NOT_MOD3:
            ++not_mod3;
            break;

          case ProteinSequenceValidity::EMPTY:
            ++empty_sequence;
            break;

          case ProteinSequenceValidity::NO_START_CODON:
            ++no_start_codon;
            break;

          case ProteinSequenceValidity::NO_STOP_CODON:
            ++no_stop_codon;
            break;

          case ProteinSequenceValidity::NONSENSE_MUTATION:
            ++nonsense_mutation;
            break;

        }


      } //for sequences

      if (not ncRNA_gene) {

        if (valid_gene) {

          ++protein_gene_count;

        } else {

          ++ill_formed_genes;

        }

      }

    } // If gene.

  }

  ExecEnv::log().info("Verified {}; Valid Protein Genes {}, Mal-formed (pseudo) {}; ncRNA Genes {}, Empty Genes: {}, Valid Sequences: {}, Not Mod3: {}, Bad Structure: {}",
                      contigId(),
                      protein_gene_count,
                      ill_formed_genes,
                      ncRNA_gene_count,
                      empty_genes,
                      valid_sequence,
                      not_mod3,
                      (no_start_codon + no_stop_codon + nonsense_mutation));

}


bool kgl::ContigReference::verifyGene(const std::shared_ptr<const GeneFeature>& gene_ptr) {

  auto coding_seq_ptr = kgl::GeneFeature::getTranscriptionSequences(gene_ptr);
  if (coding_seq_ptr->empty()) {

    return true;   // OK, loci and regions are marked as genes but have no coding sequences.

  }

  auto valid_sequence_array_ptr = coding_seq_ptr->validTranscriptionArray();

  // Valid gene should have at least 1 valid transcription.
  return not valid_sequence_array_ptr->empty();

}


bool kgl::Feature::verifyMod3(const TranscriptionFeatureMap& feature_map) const {

  bool result = true;
// Check the combined sequence length is mod 3 = 0

  ContigSize_t coding_sequence_length = 0;
  for (auto const& [offset, feature_ptr] : feature_map) {

    coding_sequence_length += feature_ptr->sequence().length();

  }

  if ((coding_sequence_length % Codon::CODON_SIZE) != 0) {

    result = false;

  }

  return result;

}

bool kgl::Feature::verifyStrand(const TranscriptionFeatureMap& Feature_map) const {

  bool result = true;

// Check the strand is consistent and not unknown.
  for (auto const& [offset, feature_ptr] : Feature_map) {

    if (feature_ptr->sequence().strand() != sequence().strand()) {

      ExecEnv::log().error("Feature: {} strand: {}, Parent id: {} strand: {} mis-match",
                           feature_ptr->id(),
                           feature_ptr->sequence().strandText(),
                           id(),
                           sequence().strandText());
      result = false;

    }

  }

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig_ref_ptr.
kgl::ProteinSequenceValidity kgl::ContigReference::checkValidProteinSequence(const AminoSequence& amino_sequence) const {

  if (not coding_table_.checkStartCodon(amino_sequence)) {

    return ProteinSequenceValidity::NO_START_CODON;

  }

  if (coding_table_.checkNonsenseMutation(amino_sequence) != 0) {

    return ProteinSequenceValidity::NONSENSE_MUTATION;

  }

  if (not coding_table_.checkStopCodon(amino_sequence)) {

    return ProteinSequenceValidity::NO_STOP_CODON;

  }

  return ProteinSequenceValidity::VALID;

}


size_t kgl::ContigReference::proteinSequenceSize(const AminoSequence& amino_sequence) const {

  switch(checkValidProteinSequence(amino_sequence)) {

    case ProteinSequenceValidity::NO_START_CODON:
      return 0;

    case ProteinSequenceValidity::NONSENSE_MUTATION:
      return coding_table_.checkNonsenseMutation(amino_sequence);

    case ProteinSequenceValidity::NO_STOP_CODON:
      return 0;

    case ProteinSequenceValidity::VALID:
      return amino_sequence.length();

    default:
      return 0;

  }

  return 0; // Never reached, to keep the compiler happy.

}

