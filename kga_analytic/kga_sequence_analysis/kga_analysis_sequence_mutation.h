//
// Created by kellerberrin on 6/07/23.
//

#ifndef KGL_DB_MUTATION_H
#define KGL_DB_MUTATION_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kgl_mutation_analysis.h"
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


using GeneContigMap = std::map<std::shared_ptr<const GeneFeature>, ContigId_t>;
class MutateGenes {

public:

  explicit MutateGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) : genome_ptr_(genome_ptr) {

    initializeGeneContigMap(genome_ptr_);

  }
  ~MutateGenes() = default;
  // Return the id of all genes from a particular contig_ref_ptr.
  std::vector<std::shared_ptr<const GeneFeature>> contigGenes(const ContigId_t& contig_id) const;
  // Mutate the population.
  void mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr);
  // Access the mutation analysis object.
  const MutateAnalysis& mutateAnalysis() const { return mutate_analysis_; }

private:

  std::shared_ptr<const GenomeReference> genome_ptr_;
  GeneContigMap gene_contig_map_;
  mutable MutateAnalysis mutate_analysis_;

  void initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr);

  // Process a gene and transcript.
  void mutateTranscript(const std::shared_ptr<const GeneFeature>& gene_ptr,
                        const FeatureIdent_t& transcript_id,
                        const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                        const std::shared_ptr<const PopulationDB>& population,
                        const std::shared_ptr<const GenomeReference>& reference_genome_ptr) const;



  // .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
  MutateStats mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                             const FeatureIdent_t& transcript_id,
                             const std::shared_ptr<const PopulationDB>& gene_population_ptr,
                             const std::shared_ptr<const GenomeReference>& reference_genome_ptr) const;

  // Multi-threaded function, statistics and other objects returned in RegionReturn.
  static std::pair<SequenceStats,bool> genomeTranscriptMutation( const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                                               const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                               const FeatureIdent_t& transcript_id,
                                                               const std::shared_ptr<const GenomeReference>& genome_ref_ptr);


};



} // Namespace.

#endif //KGL_DB_MUTATION_H
