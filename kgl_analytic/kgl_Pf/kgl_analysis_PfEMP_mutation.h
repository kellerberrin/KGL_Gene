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

private:

  std::shared_ptr<const GenomeReference> genome_ptr_;
  GeneContigMap gene_contig_map_;

  void initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr);

};


} // Namespace.

#endif //KGL_KGL_ANALYSIS_PFEMP_MUTATION_H
