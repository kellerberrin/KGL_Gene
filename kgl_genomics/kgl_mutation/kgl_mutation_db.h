//
// Created by kellerberrin on 6/07/23.
//

#ifndef KGL_DB_MUTATION_H
#define KGL_DB_MUTATION_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kgl_mutation_analysis.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Used to find multiple different variants at a particular Genome/Contig/Offset.
// Multiple variants indicate different minor alleles at the same offset.
// P.Falciparum is haploid at the blood stage, so this indicates complexity of infection
// or an issue with the generation of the VCF file.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GeneContigMap = std::map<std::shared_ptr<const GeneFeature>, ContigId_t>;
class MutateGenes {

public:

  explicit MutateGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) : genome_ptr_(genome_ptr) {

    initializeGeneContigMap(genome_ptr_);

  }
  ~MutateGenes() = default;
  // Return the id of all genes from a particular contig.
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
                        const std::shared_ptr<const GenomeReference>& reference_genome) const;

  // .first total variants, .second multiple (duplicate) variants per offset.
  static std::tuple<std::string, size_t, size_t> threadMutation(std::shared_ptr<const GenomeDB> genome_ptr,
                                                                const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                                const FeatureIdent_t& transcript_id);
  // .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
  std::tuple<size_t, size_t, size_t> mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                    const FeatureIdent_t& transcript_id,
                                                    const std::shared_ptr<const PopulationDB>& gene_population_ptr) const;

};


} // Namespace.

#endif //KGL_DB_MUTATION_H
