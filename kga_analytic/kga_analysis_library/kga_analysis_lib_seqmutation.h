//
// Created by kellerberrin on 6/07/23.
//

#ifndef KGL_DB_MUTATION_H
#define KGL_DB_MUTATION_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kga_analysis_sequence_statistics.h"
#include "kgl_mutation_variant_map.h"
#include "kgl_mutation_interval.h"
#include "kgl_mutation_transcript.h"

namespace kellerberrin::genome::analysis {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Used to find multiple different variants at a particular Genome/Contig/Offset.
// Multiple variants indicate different minor alleles at the same offset.
// P.Falciparum is haploid at the blood stage, so this indicates complexity of infection
// or an issue with the generation of the VCF file.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simple struct to return generated sequence stats.
struct SequenceStats {

  FilteredVariantStats filter_statistics_;
  CodingSequenceValidity original_sequence_{CodingSequenceValidity::VALID_PROTEIN};
  size_t original_amino_size_{0};
  CodingSequenceValidity modified_sequence_{CodingSequenceValidity::VALID_PROTEIN};
  size_t modified_amino_size_{0};

};


class MutateGenes {

public:

  MutateGenes(const std::shared_ptr<const GenomeReference>& genome_ptr, SeqVariantFilterType filtertype)
              : genome_ptr_(genome_ptr), filtertype_(filtertype) {}
  ~MutateGenes() = default;
  // Return the id of all genes from a particular contig_ref_ptr.
  [[nodiscard]] std::vector<std::shared_ptr<const GeneFeature>> contigGenes(const ContigId_t& contig_id) const;
  // Mutate the population.
  void mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr);
  // Access the mutation analysis object.
  [[nodiscard]] const MutateAnalysis& mutateAnalysis() const { return mutate_analysis_; }

private:

  std::shared_ptr<const GenomeReference> genome_ptr_;
  MutateAnalysis mutate_analysis_;
  const SeqVariantFilterType filtertype_;


  // Process a gene and transcript.
  TranscriptMutateRecord mutateTranscript(const std::shared_ptr<const GeneFeature>& gene_ptr,
                                          const FeatureIdent_t& transcript_id,
                                          const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                          const std::shared_ptr<const PopulationDB>& population,
                                          const std::shared_ptr<const GenomeReference>& reference_genome_ptr);

  // .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
  MutateStats mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                             const FeatureIdent_t& transcript_id,
                             const std::shared_ptr<const PopulationDB>& gene_population_ptr,
                             const std::shared_ptr<const GenomeReference>& reference_genome_ptr);

  // Multi-threaded function, statistics and other objects returned in RegionReturn.
  std::pair<SequenceStats,bool> genomeTranscriptMutation( const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                                          const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                          const FeatureIdent_t& transcript_id,
                                                          const std::shared_ptr<const GenomeReference>& genome_ref_ptr);


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MutateGenesReport {

public:

  MutateGenesReport(const std::shared_ptr<const GenomeReference>& genome_ref_ptr,
                    SeqVariantFilterType filter_type,
                    std::string report_directory)
  : mutate_genes_(genome_ref_ptr, filter_type), report_directory_(std::move(report_directory)) {}
  ~MutateGenesReport() = default;

  void mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr);
  void printMutateReports();

private:

  MutateGenes mutate_genes_; // Perform transcript level mutations for all genomes.
  std::string report_directory_;
  constexpr static const std::string VARIANT_COUNT_EXT_{".csv"};

};



} // Namespace.

#endif //KGL_DB_MUTATION_H
