//
// Created by kellerberrin on 28/12/20.
//

#ifndef KGL_GENOME_CONTIG_H
#define KGL_GENOME_CONTIG_H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "kgl_properties.h"
#include "kgl_genome_types.h"
#include "kgl_sequence_amino.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig_feature.h"
#include "kgl_genome_contig_aux.h"
#include "kgl_gaf_parser.h"
#include "kel_interval_unsigned.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference - A contiguous region, the associated sequence, and all features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigReference {

public:

  ContigReference(ContigId_t contig_id,
                  const std::shared_ptr<const DNA5SequenceLinear>& sequence_ptr)
                  : contig_id_(std::move(contig_id)), sequence_ptr_(sequence_ptr) {}
  ContigReference(const ContigReference&) = default;
  ~ContigReference() = default;

  ContigReference& operator=(const ContigReference&) = default;

  // Add parsed features to the different feature structures.
  [[nodiscard]] bool addContigFeature(std::shared_ptr<Feature>& feature_ptr);

  [[nodiscard]] std::vector<std::shared_ptr<const Feature>> findFeatureId(const FeatureIdent_t& feature_id) const;

  [[nodiscard]] const GeneMap& getGeneMap() const { return gene_exon_features_.geneMap(); }

  [[nodiscard]] bool setTranslationTable(const std::string& table_name) { return coding_table_.settranslationTable(table_name); }

  [[nodiscard]] std::string translationTableName() const { return coding_table_.translationTableName(); }

  void description(const std::string& desc) { description_ = desc; }

  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] const std::string& description() const { return description_; }
  [[nodiscard]] const DNA5SequenceLinear& sequence() const { return *sequence_ptr_; }
  [[nodiscard]] const TranslateToAmino& codingTable() const { return coding_table_; }
  [[nodiscard]] std::shared_ptr<const DNA5SequenceLinear> sequence_ptr() const { return sequence_ptr_; }
  [[nodiscard]] ContigSize_t contigSize() const { return sequence_ptr_->length(); }


  [[nodiscard]] static bool verifyGene(const std::shared_ptr<const GeneFeature>& gene_ptr);
  // Returns the protein sequence as the distance in amino acids between the start codon and stop codon.
  // No start and stop codon returns 0.
  [[nodiscard]] size_t proteinSequenceSize(const AminoSequence& amino_sequence) const;

  // Given a gene id and an mRNA (sequence id) return the CDS coding sequence.
  [[nodiscard]] std::optional<std::shared_ptr<const TranscriptionSequence>>
    getCodingSequence(const FeatureIdent_t& gene_id, const FeatureIdent_t& transcript_id) const;

  // Given a transcript return the associated coding sequence.
  [[nodiscard]] std::optional<DNA5SequenceCoding>
    codingSequence( const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) const;

  // Generate Amino acid sequences using the table specified for this contig_ref_ptr.
  [[nodiscard]] AminoSequence getAminoSequence(const DNA5SequenceCoding& sequence_ptr) const;
  [[nodiscard]] AminoAcid::Alphabet getAminoAcid(const Codon& codon) const { return coding_table_.getAmino(codon); }

  // Compare reference Contigs- mainly used for testing.
  [[nodiscard]] bool equivalent(const ContigReference& compare_contig) const;
  // Check start codon, stop codon codon and nonsense mutation.
  [[nodiscard]] ProteinSequenceValidity checkValidProteinSequence(const AminoSequence& amino_sequence) const;
  // Check all of the above and Mod3.
  [[nodiscard]] ProteinSequenceValidity checkValidCodingSequence(const DNA5SequenceCoding& coding_sequence) const;

  // Wire-up the contig_ref_ptr features
  void verifyFeatureHierarchy();
  void verifyCDSPhasePeptide();

private:

  ContigId_t contig_id_;
  std::string description_;
  std::shared_ptr<const DNA5SequenceLinear> sequence_ptr_;  // The contig_ref_ptr unstranded DNA sequence.
  GeneExonFeatures gene_exon_features_;  // All the genes and sequences defined for this contig_ref_ptr.
  TranslateToAmino coding_table_;  // Amino Acid translation table, unique for contig_ref_ptr (e.g. mitochondria)

};




}   // end namespace



#endif //KGL_GENOME_CONTIG_H
