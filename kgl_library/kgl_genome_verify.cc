//
// Created by kellerberrin on 12/11/17.
//


#include "kel_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void kgl::ContigFeatures::verifyCDSPhasePeptide() {

  // Iterate through all the features looking for Genes.
  size_t gene_count = 0;
  size_t well_formed_sequences = 0;
  size_t ill_formed_sequences = 0;
  size_t empty_genes = 0;

  ExecEnv::log().info("ContigFeatures::verifyCDSPhasePeptide; Verifying {} Genes using amino translation table: {}", contigId(), coding_table_.translationTableName());

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

          ExecEnv::log().warn("Gene : {} Offset: {} Problem verifying CDS structure - gene sub-features print out",
                              gene_ptr->id(),
                              gene_ptr->sequence().begin());
          feature.second->recusivelyPrintsubfeatures();

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


bool kgl::ContigFeatures::verifyCodingSequences(const std::shared_ptr<const GeneFeature> gene_ptr,
                                                const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const {

  bool result = true;

  if (coding_seq_ptr->size() == 0) {

    ExecEnv::log().error("codingSequence(), empty CodingSequenceArray");

  }

  for (const auto& sequence : coding_seq_ptr->getMap()) {

    if (sequence.second->getSortedCDS().empty()) {

      ExecEnv::log().error("codingSequence(), no corresponding CDS features found");
      continue;

    }

    std::shared_ptr<DNA5SequenceCoding> coding_sequence_ptr = sequence_ptr_->codingSequence(sequence.second);

    if (not coding_table_.checkStartCodon(coding_sequence_ptr)) {

      std::vector<std::string> description_vec;
      if (not gene_ptr->getAttributes().getDescription(description_vec)) {

        ExecEnv::log().error("Cannot get description vector for Gene: {}", gene_ptr->id());
        description_vec.clear();

      }

      if (description_vec.empty()) {

        description_vec.emplace_back("No Description");

      }
      std::string gene_description = description_vec.front();

      ExecEnv::log().vwarn("No START codon Gene: {}: {}, Sequence (mRNA): {} | first codon: {}",
                           sequence.second->getGene()->id(),
                           gene_description,
                           sequence.second->getCDSParent()->id(),
                           coding_table_.firstCodon(coding_sequence_ptr).getSequenceAsString());
//      gene_ptr->recusivelyPrintsubfeatures();
      result = false;
    }
    if (not coding_table_.checkStopCodon(coding_sequence_ptr)) {

      std::vector<std::string> description_vec;

      if (not gene_ptr->getAttributes().getDescription(description_vec)) {

        ExecEnv::log().error("Cannot get description vector for Gene: {}", gene_ptr->id());
        description_vec.clear();

      }

      if (description_vec.empty()) {

        description_vec.emplace_back("No Description");

      }
      std::string gene_description = description_vec.front();

      ExecEnv::log().vwarn("No STOP codon: {} Gene: {}: {}, Sequence (mRNA): {} | last codon: {}",
                           (Codon::codonLength(coding_sequence_ptr)-1),
                           sequence.second->getGene()->id(),
                           gene_description,
                           sequence.second->getCDSParent()->id(),
                           coding_table_.lastCodon(coding_sequence_ptr).getSequenceAsString());
//      gene_ptr->recusivelyPrintsubfeatures();
      result = false;
    }
    size_t nonsense_index = coding_table_.checkNonsenseMutation(coding_sequence_ptr);
    if (nonsense_index > 0) {

      std::vector<std::string> description_vec;
      if (not gene_ptr->getAttributes().getDescription(description_vec)) {

        ExecEnv::log().error("Cannot get description vector for Gene: {}", gene_ptr->id());
        description_vec.clear();

      }

      if (description_vec.empty()) {

        description_vec.emplace_back("No Description");

      }
      std::string gene_description = description_vec.front();

      ExecEnv::log().vwarn("NONSENSE mutation codon:{} Gene: {}: {}, Sequence (mRNA): {} | stop codon: {}",
                           nonsense_index,
                           sequence.second->getGene()->id(),
                           gene_description,
                           sequence.second->getCDSParent()->id(),
                           Codon(coding_sequence_ptr, nonsense_index).getSequenceAsString());
//      gene_ptr->recusivelyPrintsubfeatures();
      result = false;
    }

  } // for cds group

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig.
bool kgl::ContigFeatures::verifyDNACodingSequence(std::shared_ptr<const DNA5SequenceCoding> coding_sequence_ptr) const {

  bool result = coding_table_.checkStopCodon(coding_sequence_ptr);
  result = result and coding_table_.checkStartCodon(coding_sequence_ptr);
  result = result and coding_table_.checkNonsenseMutation(coding_sequence_ptr) == 0;

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig.
bool kgl::ContigFeatures::verifyProteinSequence(std::shared_ptr<const AminoSequence> amino_sequence_ptr) const {

  bool result = coding_table_.checkStopCodon(amino_sequence_ptr);
  result = result and coding_table_.checkStartCodon(amino_sequence_ptr);
  result = result and coding_table_.checkNonsenseMutation(amino_sequence_ptr) == 0;

  return result;

}


// Verifies a coding sequence using the amino coding table defined for the contig.
kgl::ProteinSequenceAnalysis kgl::ContigFeatures::proteinSequenceAnalysis(std::shared_ptr<const AminoSequence> amino_sequence_ptr) const {

  if (not coding_table_.checkStartCodon(amino_sequence_ptr)) {

    return ProteinSequenceAnalysis::NO_START_CODON;

  }

  if (coding_table_.checkNonsenseMutation(amino_sequence_ptr) != 0) {

    return ProteinSequenceAnalysis::NONSENSE_MUTATION;

  }

  if (not coding_table_.checkStopCodon(amino_sequence_ptr)) {

    return ProteinSequenceAnalysis::NO_STOP_CODON;

  }

  return ProteinSequenceAnalysis::VALID_SEQUENCE;

}


size_t kgl::ContigFeatures::proteinSequenceSize(std::shared_ptr<const AminoSequence> amino_sequence_ptr) const {

  switch(proteinSequenceAnalysis(amino_sequence_ptr)) {

    case ProteinSequenceAnalysis::NO_START_CODON:
      return 0;

    case ProteinSequenceAnalysis::NONSENSE_MUTATION:
      return coding_table_.checkNonsenseMutation(amino_sequence_ptr);

    case ProteinSequenceAnalysis::NO_STOP_CODON:
      return 0;

    case ProteinSequenceAnalysis::VALID_SEQUENCE:
      return amino_sequence_ptr->length();

  }

  return 0; // Never reached, to keep the compiler happy.

}


