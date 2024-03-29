//
// Created by kellerberrin on 17/11/17.
//

#ifndef KGL_PHYLOGENETIC_ANALYSIS_H
#define KGL_PHYLOGENETIC_ANALYSIS_H


#include <memory>
#include <fstream>
#include "kgl_variant_db_population.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_io_gff_fasta.h"
#include "kgl_pfgenome_aux.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_sequence_compare_impl.h"




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object implements a series of high level functions for the kgl_classification application.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization level namespace



class MutationItem {

public:

  MutationItem() = default;
  MutationItem(const MutationItem&) = default;
  ~MutationItem() = default;

  MutationItem& operator=(const MutationItem&) = default;

  EditItem DNA_mutation;
  std::string reference_codon;
  std::string mutation_codon;
  EditItem amino_mutation;
  ContigId_t contig_id;
  ContigOffset_t contig_offset;

  [[nodiscard]] std::string mapKey() const {

    std::stringstream ss;
    ss << amino_mutation.reference_offset << amino_mutation.reference_char << amino_mutation.mutant_char;
    ss << DNA_mutation.reference_offset << DNA_mutation.reference_char << DNA_mutation.mutant_char;

    return ss.str(); }

};


class MutationEditVector {

public:

  MutationEditVector() = default;
  MutationEditVector(const MutationEditVector&) = default;
  ~MutationEditVector() = default;

  std::vector<MutationItem> mutation_vector;

  [[nodiscard]]  bool hasIndel() const {

    for (auto mutation : mutation_vector) {

      if (mutation.DNA_mutation.isIndel()) return true;

    }

    return false;

  }

};


enum class SequenceAnalysisType { DNA, VARIANT, SNP, SIZE, ENTROPY, LEMPEL_ZIV};

class GenomicMutation {

public:

  GenomicMutation() = default;
  ~GenomicMutation() = default;

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  [[nodiscard]]  static bool writeMutantProteins( const std::string& fastaFile,
                                                  const std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>>& amino_seq_vector);

  [[nodiscard]]  static bool writeMutantDNA( const std::string& fasta_file,
                                             const std::vector<std::pair<std::string, std::shared_ptr<DNA5SequenceLinear>>>& dna_seq_vector);

  // If sequence name = "" then the code reads all sequences in the fasta file.
  [[nodiscard]]  static bool readFastaProteins( const std::string& fasta_file,
                                                std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>>& amino_seq_vector);

  [[nodiscard]] static std::optional<std::pair<DNA5SequenceCoding, DNA5SequenceCoding>>
                compare5Prime( const ContigId_t& contig_id,
                               const FeatureIdent_t& gene_id,
                               const FeatureIdent_t& transcript_id,
                               ContigSize_t region_size,
                               const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                               const std::shared_ptr<const GenomeDB>& genome_db_ptr);

  [[nodiscard]] static std::optional<std::pair<DNA5SequenceCoding, DNA5SequenceCoding>>
                compare3Prime( const ContigId_t& contig_id,
                               const FeatureIdent_t& gene_id,
                               const FeatureIdent_t& transcript_id,
                               ContigSize_t region_size,
                               const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                               const std::shared_ptr<const GenomeDB>& genome_db_ptr);


  [[nodiscard]] static bool outputRegionCSV( const std::string &file_name,
                                             LinearDistanceMetric linear_dna_distance_metric,
                                             CodingDistanceMetric coding_dna_distance_metric,
                                             AminoDistanceMetric amino_distance_metric,
                                             const std::shared_ptr<const GenomeReference>& genome_db,
                                             const std::shared_ptr<const PopulationDB>& pop_variant_ptr);


  [[nodiscard]]  static bool outputDNASequenceCSV( const std::string &file_name,
                                                   SequenceAnalysisType analysis_type,
                                                   CodingDistanceMetric dna_distance_metric,
                                                   const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                                   const std::shared_ptr<const PopulationDB>& pop_variant_ptr);

  [[nodiscard]]  static bool outputAminoSequenceCSV( const std::string &file_name,
                                                     AminoDistanceMetric amino_distance_metric,
                                                     const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                                     const std::shared_ptr<const PopulationDB>& pop_variant_ptr);


  [[nodiscard]]  static bool outputAminoMutationCSV( const std::string &file_name,
                                                     const ContigId_t& contig_id,
                                                     const FeatureIdent_t& gene_id,
                                                     const FeatureIdent_t& transcript_id,
                                                     const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                                     const std::shared_ptr<const PopulationDB>& pop_variant_ptr);

  [[nodiscard]]  static bool outputDNAMutationCSV( const std::string &file_name,
                                                   const ContigId_t& contig_id,
                                                   const FeatureIdent_t& gene_id,
                                                   const FeatureIdent_t& transcript_id,
                                                   const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                                   const std::shared_ptr<const PopulationDB>& pop_variant_ptr,
                                                   const PfGenomeAuxData& aux_Pf3k_data);

private:

  [[nodiscard]] static std::string outputSequenceHeader(char delimiter, const std::shared_ptr<const PopulationDB>& pop_variant_ptr);
  [[nodiscard]] static std::string outputRegionHeader(char delimiter);
  [[nodiscard]] static std::string outputSequence( char delimiter,
                                                   CodingDistanceMetric dna_distance_metric,
                                                   AminoDistanceMetric amino_distance_metric,
                                                   const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                                   const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                                                   const std::shared_ptr<const GenomeDB>& genome_db_ptr);


};



}   // end namespace genome




#endif //KGL_PHYLOGENETIC_ANALYSIS_H
