//
// Created by kellerberrin on 6/07/23.
//

#ifndef KGL_ANALYSIS_PFEMP_MUTATION_H
#define KGL_ANALYSIS_PFEMP_MUTATION_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace


using GeneContigMap = std::map<std::shared_ptr<const GeneFeature>, ContigId_t>;
class MutateGenes {

public:

  explicit MutateGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) : genome_ptr_(genome_ptr) {

    initializeGeneContigMap(genome_ptr_);

  }
  ~MutateGenes() = default;
  // Return the id of all genes from a particular contig.
  std::vector<std::shared_ptr<const GeneFeature>> contigGenes(const ContigId_t& contig_id) const;
  // Process a gene and transcript.
  void mutateTranscript(const std::shared_ptr<const GeneFeature>& gene_ptr,
                        const FeatureIdent_t& transcript,
                        const std::shared_ptr<const PopulationDB>& population,
                        const std::shared_ptr<const GenomeReference>& reference_genome) const;



private:

  std::shared_ptr<const GenomeReference> genome_ptr_;
  GeneContigMap gene_contig_map_;

  void initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr);

  static bool threadMutation( std::shared_ptr<const GenomeDB> genome_ptr,
                              const std::shared_ptr<const GeneFeature>& gene_ptr,
                              const FeatureIdent_t& transcript_id);
  void mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                      const FeatureIdent_t& transcript_id,
                      const std::shared_ptr<const PopulationDB>& gene_population_ptr) const;

};


} // Namespace.

#endif //KGL_KGL_ANALYSIS_PFEMP_MUTATION_H
