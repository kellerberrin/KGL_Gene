//
// Created by kellerberrin on 20/6/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
#define KGL_ANALYSIS_MUTATION_GENE_ALLELE_H


#include "kgl_analysis_mutation_gene.h"



namespace kellerberrin::genome {   //  organization::project level namespace



class GenerateGeneAllele {

public:

  explicit GenerateGeneAllele(const std::shared_ptr<const PopulationDB>& unphased_population_ptr)
  : unphased_population_ptr_(unphased_population_ptr) {

    sorted_allele_map_ = GenomeMutation::ensemblIndex(unphased_population_ptr_);

  }
  ~GenerateGeneAllele() = default;

  void generateVariantList(const std::vector<std::string>& gene_list) const;

private:

  std::shared_ptr<const PopulationDB> unphased_population_ptr_;
  std::shared_ptr<const EnsemblIndexMap> sorted_allele_map_;

};



} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
