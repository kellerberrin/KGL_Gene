//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_genome_feature.h"
#include "kgl_genome_contig.h"
#include "kgl_genome_verify.h"


namespace kgl = kellerberrin::genome;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple object to hold coding sequence validity statistics.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::SequenceValidityStatistics::updateValidity(CodingSequenceValidity validity, size_t amino_size) {

  switch(validity) {

    case CodingSequenceValidity::NCRNA:
      ++ncRNA_;
      break;

    case CodingSequenceValidity::VALID_PROTEIN:
      ++valid_protein_;
      break;

    case CodingSequenceValidity::EMPTY:
      ++empty_;
      break;

    case CodingSequenceValidity::NOT_MOD3:
      ++not_mod3_;
      break;

    case CodingSequenceValidity::NO_START_CODON:
      ++no_start_codon_;
      break;

    case CodingSequenceValidity::NONSENSE_MUTATION:
      ++nonsense_mutation_;
      nonsense_mutation_size_.push_back(amino_size);
      break;

    case CodingSequenceValidity::NO_STOP_CODON:
      ++no_stop_codon_;
      break;

  }

}


void kgl::SequenceValidityStatistics::updateTranscriptArray(const std::shared_ptr<const TranscriptionSequenceArray>& transcript_array_ptr) {

  for (auto const& [trnascript_id, transcript_ptr]: transcript_array_ptr->getMap()) {

    updateTranscript(transcript_ptr);

  }

}


void kgl::SequenceValidityStatistics::updateTranscript(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) {

  auto sequence_validity = TranscriptionSequence::checkSequenceStatus(transcript_ptr);
  updateValidity(sequence_validity)  ;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::ContigReference::verifyGeneFeatures() {

  // Iterate through all the features looking for Genes.
  SequenceValidityStatistics validity_statistics;

  for(const auto& [feature_id, feature_ptr  ] : gene_exon_features_.offsetFeatureMap()) {

    if (feature_ptr->isGene()) {

      auto gene_ptr = std::static_pointer_cast<const GeneFeature>(feature_ptr);
      auto sequence_array_ptr = kgl::GeneFeature::getTranscriptionSequences(gene_ptr);

      validity_statistics.updateTranscriptArray(sequence_array_ptr);


    } // If gene.

  } // For all features.

  ExecEnv::log().info("Verified {}; Sequences: {}, ncRNA {}, Valid Protein {}, Mal-formed (pseudo) {};  Empty Genes: {}, Not Mod3: {}, Bad Structure: {}",
                      contigId(),
                      validity_statistics.totalSequence(),
                      validity_statistics.ncRNA(),
                      validity_statistics.validProtein(),
                      validity_statistics.invalidProtein(),
                      validity_statistics.empty(),
                      validity_statistics.notMod3(),
                      (validity_statistics.invalidProtein() - validity_statistics.notMod3()));


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

  if (Codon::codonRemainder(coding_sequence_length) != 0) {

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
kgl::CodingSequenceValidity kgl::ContigReference::checkValidProteinSequence(const AminoSequence& amino_sequence) const {

  if (not coding_table_.checkStartCodon(amino_sequence)) {

    return CodingSequenceValidity::NO_START_CODON;

  }

  auto [first_size, first_found] = coding_table_.firstStopSequenceSize(amino_sequence);
  if (first_size != amino_sequence.length()) {

    return CodingSequenceValidity::NONSENSE_MUTATION;

  }

  if (not coding_table_.checkStopCodon(amino_sequence)) {

    return CodingSequenceValidity::NO_STOP_CODON;

  }

  return CodingSequenceValidity::VALID_PROTEIN;

}


kgl::CodingSequenceValidity
kgl::ContigReference::checkValidCodingSequence(const DNA5SequenceCoding& coding_sequence) const {


  auto sequence_length = coding_sequence.length();
  if (Codon::codonRemainder(sequence_length) != 0) {

    return CodingSequenceValidity::NOT_MOD3;

  }

  auto amino_sequence = getAminoSequence(coding_sequence);
  return checkValidProteinSequence(amino_sequence);

}



std::pair<kgl::CodingSequenceValidity, size_t> kgl::ContigReference::proteinSequenceSize(const AminoSequence& amino_sequence) const {

  auto validity_code = checkValidProteinSequence(amino_sequence);
  auto [first_size, first_found] = coding_table_.firstStopSequenceSize(amino_sequence);

  return {validity_code, first_size}; // Never reached, to keep the compiler happy.

}


std::pair<kgl::CodingSequenceValidity, size_t>
kgl::ContigReference::codingProteinSequenceSize(const DNA5SequenceCoding& coding_sequence) const {


  auto sequence_length = coding_sequence.length();
  if (Codon::codonRemainder(sequence_length) != 0) {

    return { CodingSequenceValidity::NOT_MOD3, Codon::codonLength(coding_sequence)};

  }

  auto amino_sequence = getAminoSequence(coding_sequence);
  return proteinSequenceSize(amino_sequence);

}

