//
// Created by kellerberrin on 17/11/17.
//

#ifndef KGL_PHYLOGENETIC_ANALYSIS_H
#define KGL_PHYLOGENETIC_ANALYSIS_H


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_db.h"
#include "kgl_filter.h"
#include "kgl_gff_fasta.h"
#include "kgl_genome_aux_csv.h"
#include "kgl_sequence_compare_impl.h"
#include "kgl_variant_phasing_statistics.h"




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object implements a series of high level functions for the kgl_phylogenetic application.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




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

  std::string mapKey() const {

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

  bool hasIndel() const {

    for (auto mutation : mutation_vector) {

      if (mutation.DNA_mutation.isIndel()) return true;

    }

    return false;

  }

};


class ApplicationAnalysis {

public:

  ApplicationAnalysis() = default;
  virtual ~ApplicationAnalysis() = default;

  // Sequences are presented as a pair of a sequence name and an amino sequences.
  static bool writeMutantProteins(const std::string& fastaFile,
                                  const std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>>& amino_seq_vector);

  static bool writeMutantDNA(const std::string& fasta_file,
                             const std::vector<std::pair<std::string, std::shared_ptr<DNA5SequenceLinear>>>& dna_seq_vector);

  // If sequence name = "" then the code reads all sequences in the fasta file.
  static bool readFastaProteins(const std::string& fasta_file,
                                std::vector<std::pair<std::string, std::shared_ptr<AminoSequence>>>& amino_seq_vector);

  static bool compare5Prime(const ContigId_t& contig_id,
                            const FeatureIdent_t& gene_id,
                            const FeatureIdent_t& sequence_id,
                            ContigSize_t region_size,
                            const std::shared_ptr<const GenomeDatabase>& genome_db,
                            const std::shared_ptr<const GenomeVariant>& genome_variant,
                            std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                            std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector);

  static bool compare3Prime(const ContigId_t& contig_id,
                            const FeatureIdent_t& gene_id,
                            const FeatureIdent_t& sequence_id,
                            ContigSize_t region_size,
                            const std::shared_ptr<const GenomeDatabase>& genome_db,
                            const std::shared_ptr<const GenomeVariant>& genome_variant,
                            std::shared_ptr<DNA5SequenceCoding>& reference_sequence,
                            std::vector<std::shared_ptr<DNA5SequenceCoding>>& mutant_sequence_vector);

  static bool outputSequenceCSV(const std::string &file_name,
                                std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric,
                                std::shared_ptr<const GenomeDatabase> genome_db,
                                std::shared_ptr<const PhasedPopulation> pop_variant_ptr);

  static bool outputAminoMutationCSV(const std::string &file_name,
                                     const ContigId_t& contig_id,
                                     const FeatureIdent_t& gene_id,
                                     const FeatureIdent_t& sequence_id,
                                     std::shared_ptr<const GenomeDatabase> genome_db,
                                     std::shared_ptr<const PhasedPopulation> pop_variant_ptr);

  static bool outputDNAMutationCSV(const std::string &file_name,
                                   const ContigId_t& contig_id,
                                   const FeatureIdent_t& gene_id,
                                   const FeatureIdent_t& sequence_id,
                                   std::shared_ptr<const GenomeDatabase> genome_db,
                                   std::shared_ptr<const PhasedPopulation> pop_variant_ptr,
                                   const GenomeAuxData& aux_Pf3k_data,
                                   std::shared_ptr<const PopulationPhasingStatistics> phasing_stats);


private:

  static std::string outputSequenceHeader(char delimiter);
  static std::string outputSequence(char delimiter,
                                    std::shared_ptr<const GlobalDNASequenceDistance> dna_distance_metric,
                                    std::shared_ptr<const GlobalAminoSequenceDistance> amino_distance_metric,
                                    std::shared_ptr<const CodingSequence> coding_sequence,
                                    std::shared_ptr<const GenomeDatabase> genome_db,
                                    std::shared_ptr<const GenomeVariant> genome_variant);


};



}   // namespace genome
}   // namespace kellerberrin





#endif //KGL_PHYLOGENETIC_ANALYSIS_H
