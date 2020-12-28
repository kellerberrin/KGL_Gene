//
// Created by kellerberrin on 28/12/20.
//

#ifndef KGL_GENOME_CONTIG_H
#define KGL_GENOME_CONTIG_H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <kgl_properties.h>
#include "kgl_genome_types.h"
#include "kgl_sequence_amino.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig_feature.h"
#include "kgl_genome_contig_aux.h"
#include "kgl_gaf_parser.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigReference - A contiguous region, the associated sequence, and all features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class ProteinSequenceAnalysis { VALID_SEQUENCE, NO_START_CODON, NONSENSE_MUTATION, NO_STOP_CODON };


class ContigReference {

public:

  ContigReference(const ContigId_t& contig_id,
                  const std::shared_ptr<const DNA5SequenceContig>& sequence_ptr) : contig_id_(contig_id), sequence_ptr_(sequence_ptr) {}
  ContigReference(const ContigReference&) = default;
  ~ContigReference() = default;

  ContigReference& operator=(const ContigReference&) = default;

  // Add parsed features to the different feature structures.
  [[nodiscard]] bool addGeneExonFeature(std::shared_ptr<Feature>& feature_ptr);
  [[nodiscard]] bool addAuxFeature(std::shared_ptr<Feature>& feature_ptr);

  // false if not found.
  [[nodiscard]] bool findFeatureId( const FeatureIdent_t& feature_id,
                                    std::vector<std::shared_ptr<const Feature>>& feature_ptr_vec) const {
    return gene_exon_features_.findFeatureId(feature_id, feature_ptr_vec);
  }
  // false if offset is not in a gene, else (true) returns a vector of ptrs to the genes.
  //  [[nodiscard]] bool findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const;

  [[nodiscard]] const GeneMap& getGeneMap() const { return gene_exon_features_.geneMap(); }

  // Return all Aux genome features in this contig.
  [[nodiscard]] const AuxContigFeatures& getAuxContigFeatures() const { return aux_contig_features_; }

  [[nodiscard]] bool setTranslationTable(const std::string& table_name) { return coding_table_.settranslationTable(table_name); }

  [[nodiscard]] std::string translationTableName() const { return coding_table_.translationTableName(); }

  void description(const std::string& desc) { description_ = desc; }

  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] const std::string& description() const { return description_; }
  [[nodiscard]] const DNA5SequenceContig& sequence() const { return *sequence_ptr_; }
  [[nodiscard]] std::shared_ptr<const DNA5SequenceContig> sequence_ptr() const { return sequence_ptr_; }
  [[nodiscard]] ContigSize_t contigSize() const { return sequence_ptr_->length(); }


  [[nodiscard]] bool verifyDNACodingSequence(const DNA5SequenceCoding& coding_dna) const;
  [[nodiscard]] bool verifyProteinSequence(const AminoSequence& amino_sequence) const;
  [[nodiscard]] ProteinSequenceAnalysis proteinSequenceAnalysis(const AminoSequence& amino_sequence) const;
  // Returns the protein sequence as the distance in amino acids between the start codon and stop codon.
  // No start and stop codon returns 0.
  [[nodiscard]] size_t proteinSequenceSize(const AminoSequence& amino_sequence) const;

  // Given a gene id and an mRNA (sequence id) return the CDS coding sequence.
  [[nodiscard]] bool getCodingSequence( const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const;

  // Given a CDS coding sequence, return the corresponding DNA base sequence (strand adjusted).
  [[nodiscard]] bool getDNA5SequenceCoding( const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                                            DNA5SequenceCoding& dna_coding) const;

  // Generate Amino acid sequences using the table specified for this contig.
  [[nodiscard]] AminoSequence getAminoSequence(const DNA5SequenceCoding& sequence_ptr) const;
  [[nodiscard]] AminoAcid::Alphabet getAminoAcid(const Codon& codon) const { return coding_table_.getAmino(codon); }

  // Wire-up the contig features
  void verifyFeatureHierarchy();
  void verifyAuxillaryHierarchy();
  void verifyCDSPhasePeptide();

private:

  ContigId_t contig_id_;
  std::string description_;
  std::shared_ptr<const DNA5SequenceContig> sequence_ptr_;  // The contig unstranded DNA sequence.
  GeneExonFeatures gene_exon_features_;  // All the genes and exons defined for this contig.
  AuxContigFeatures aux_contig_features_;
  TranslateToAmino coding_table_;  // Amino Acid translation table, unique for contig (e.g. mitochondria)

  // Check all gene coding sequences for start and end codons and nonsense (intermediate stop codon) mutations.
  [[nodiscard]] bool verifyCodingSequences( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                            const std::shared_ptr<const CodingSequenceArray>& coding_seq_ptr) const;

};




}   // end namespace



#endif //KGL_GENOME_CONTIG_H
