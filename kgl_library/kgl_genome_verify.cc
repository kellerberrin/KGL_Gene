//
// Created by kellerberrin on 12/11/17.
//


#include "kel_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::ContigReference::verifyCDSPhasePeptide() {

  // Iterate through all the features looking for Genes.
  size_t gene_count = 0;
  size_t well_formed_sequences = 0;
  size_t ill_formed_sequences = 0;
  size_t empty_genes = 0;

  ExecEnv::log().info("ContigReference::verifyCDSPhasePeptide; Verifying {} Genes using amino translation table: {}", contigId(), coding_table_.translationTableName());

  for(const auto& feature : gene_exon_features_.offsetFeatureMap()) {

    if(feature.second->isGene()) {

      ++gene_count;
      const std::shared_ptr<const GeneFeature> gene_ptr = std::static_pointer_cast<const GeneFeature>(feature.second);
      const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr);
      if (coding_seq_ptr->size() == 0) { // No CDS coding sequence available.

        ++empty_genes;

      }

      if (coding_seq_ptr->size() > 0) {

        if (not gene_ptr->verifyCDSPhase(coding_seq_ptr)) {

          ExecEnv::log().vinfo("Gene : {} Offset: {} Problem verifying CDS structure - gene sub-features print out",
                              gene_ptr->id(),
                              gene_ptr->sequence().begin());

          if (ExecEnv::log().verbose()) {

            feature.second->recusivelyPrintsubfeatures();

          }


          ill_formed_sequences += coding_seq_ptr->size();

        }

        if (verifyCodingSequences(gene_ptr, coding_seq_ptr)) {

          well_formed_sequences += coding_seq_ptr->size();

        } else {

          ill_formed_sequences += coding_seq_ptr->size();

        }

      }

    }

  }

  ExecEnv::log().info("Verified {} found: {} Genes with {} well-formed and {} mal-formed (pseudo) coding sequences",
                      contigId(),
                      gene_count,
                      well_formed_sequences,
                      ill_formed_sequences);

}


bool kgl::ContigReference::verifyCodingSequences(const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                 const std::shared_ptr<const CodingSequenceArray>& coding_seq_ptr) const {

  bool result = true;

  if (coding_seq_ptr->empty()) {

    ExecEnv::log().error("gene: {}; verifyCodingSequences(), empty CodingSequenceArray", gene_ptr->id());

  }

  for (const auto& sequence : coding_seq_ptr->getMap()) {

    if (sequence.second->getSortedCDS().empty()) {

      ExecEnv::log().error("gene: {}; verifyCodingSequences(), no corresponding CDS features found", gene_ptr->id());
      continue;

    }

    DNA5SequenceCoding coding_sequence = sequence_ptr_->codingSequence(sequence.second);

    if (not coding_table_.checkStartCodon(coding_sequence)) {

      ExecEnv::log().vinfo("No START codon Gene: {}, Sequence (mRNA): {} | first codon: {}",
                           sequence.second->getGene()->id(),
                           sequence.second->getCDSParent()->id(),
                           coding_table_.firstCodon(coding_sequence).getSequenceAsString());

      if (ExecEnv::log().verbose()) {

        gene_ptr->recusivelyPrintsubfeatures();

      }

      result = false;
    }
    if (not coding_table_.checkStopCodon(coding_sequence)) {


      ExecEnv::log().vinfo("No STOP codon: {} Gene: {}, Sequence (mRNA): {} | last codon: {}",
                           (Codon::codonLength(coding_sequence)-1),
                           sequence.second->getGene()->id(),
                           sequence.second->getCDSParent()->id(),
                           coding_table_.lastCodon(coding_sequence).getSequenceAsString());

      if (ExecEnv::log().verbose()) {

        gene_ptr->recusivelyPrintsubfeatures();

      }

      result = false;
    }
    size_t nonsense_index = coding_table_.checkNonsenseMutation(coding_sequence);
    if (nonsense_index > 0) {


      ExecEnv::log().vinfo("NONSENSE mutation codon:{} Gene: {}, Sequence (mRNA): {} | stop codon: {}",
                           nonsense_index,
                           sequence.second->getGene()->id(),
                           sequence.second->getCDSParent()->id(),
                           Codon(coding_sequence, nonsense_index).getSequenceAsString());
      if (ExecEnv::log().verbose()) {

        gene_ptr->recusivelyPrintsubfeatures();

      }

      result = false;

    }

  } // for cds group

  return result;

}


bool kgl::Feature::verifyCDSPhase(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const {

  bool result = true;
  // Check for mod3
  for(const auto& sorted_cds : coding_seq_ptr->getMap()) {

    result = result and verifyMod3(sorted_cds.second->getSortedCDS());
    result = result and verifyStrand(sorted_cds.second->getSortedCDS());

  }

  return result;

}


bool kgl::Feature::verifyMod3(const SortedCDS& sorted_cds) const {

  bool result = true;
// Check the combined sequence length is mod 3 = 0

  ContigSize_t coding_sequence_length = 0;
  for (auto cds : sorted_cds) {

    coding_sequence_length += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  if ((coding_sequence_length % Codon::CODON_SIZE) != 0) {

    ExecEnv::log().vwarn("Gene: {} offset: {} CDS coding sequence length mod 3 not zero : {}",
                        id(),
                        sequence().begin(),
                        (coding_sequence_length % 3));

    result = false;

  }

  return result;

}

bool kgl::Feature::verifyStrand(const SortedCDS& sorted_cds) const {

  bool result = true;

// Check the strand is consistent and not unknown.
  for (auto cds : sorted_cds) {

    if (cds.second->sequence().strand() != sequence().strand()) {

      ExecEnv::log().error("CDS: {} offset: {} strand: {}, parent sequence strand: {} mis-match",
                           cds.second->id(),
                           cds.second->sequence().begin(),
                           static_cast<char>(cds.second->sequence().strand()),
                           static_cast<char>(sequence().strand()));
      result = false;

    }

  }

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig.
bool kgl::ContigReference::verifyDNACodingSequence(const DNA5SequenceCoding& coding_sequence_ptr) const {

  bool result = coding_table_.checkStopCodon(coding_sequence_ptr);
  result = result and coding_table_.checkStartCodon(coding_sequence_ptr);
  result = result and coding_table_.checkNonsenseMutation(coding_sequence_ptr) == 0;

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig.
bool kgl::ContigReference::verifyProteinSequence(const AminoSequence& amino_sequence) const {

  bool result = coding_table_.checkStopCodon(amino_sequence);
  result = result and coding_table_.checkStartCodon(amino_sequence);
  result = result and coding_table_.checkNonsenseMutation(amino_sequence) == 0;

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig.
kgl::ProteinSequenceAnalysis kgl::ContigReference::proteinSequenceAnalysis(const AminoSequence& amino_sequence) const {

  if (not coding_table_.checkStartCodon(amino_sequence)) {

    return ProteinSequenceAnalysis::NO_START_CODON;

  }

  if (coding_table_.checkNonsenseMutation(amino_sequence) != 0) {

    return ProteinSequenceAnalysis::NONSENSE_MUTATION;

  }

  if (not coding_table_.checkStopCodon(amino_sequence)) {

    return ProteinSequenceAnalysis::NO_STOP_CODON;

  }

  return ProteinSequenceAnalysis::VALID_SEQUENCE;

}


size_t kgl::ContigReference::proteinSequenceSize(const AminoSequence& amino_sequence) const {

  switch(proteinSequenceAnalysis(amino_sequence)) {

    case ProteinSequenceAnalysis::NO_START_CODON:
      return 0;

    case ProteinSequenceAnalysis::NONSENSE_MUTATION:
      return coding_table_.checkNonsenseMutation(amino_sequence);

    case ProteinSequenceAnalysis::NO_STOP_CODON:
      return 0;

    case ProteinSequenceAnalysis::VALID_SEQUENCE:
      return amino_sequence.length();

  }

  return 0; // Never reached, to keep the compiler happy.

}


