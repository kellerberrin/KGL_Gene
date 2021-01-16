//
// Created by kellerberrin on 15/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_H
#define KGL_ANALYSIS_MUTATION_GENE_H

#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GeneMutation {

  GenomeId_t genome;
  ContigId_t contig;
  FeatureIdent_t gene_id;
  std::string gene_name;
  std::string description;
  std::string biotype;
  bool valid_protein{false};
  std::string gaf_id;
  ContigOffset_t gene_begin{0};
  ContigOffset_t gene_end{0};
  ContigSize_t gene_size{0};
  std::string strand;
  size_t exons{0};
  size_t attribute_size{0};
  size_t variant_count{0};  // Total number of genomes found
  size_t genome_count{0};   // Total number of genomes.
  size_t genome_variant{0};  // Number of genomes that contain variants for this gene.
  size_t homozygous{0};
  size_t heterozygous{0};

};


class GenomeMutation {

public:

  GenomeMutation() = default;
  ~GenomeMutation() = default;

  // This analysis is performed first
  bool genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_reference);
  // Then this analysis.
  bool variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);

  [[nodiscard]] const std::vector<GeneMutation>& geneVector() const { return gene_vector_; }
  [[nodiscard]] std::vector<GeneMutation>& geneVector() { return gene_vector_; }


private:

  std::vector<GeneMutation> gene_vector_;

};



} // namespace




#endif //KGL_KGL_ANALYSIS_MUTATION_GENE_H
