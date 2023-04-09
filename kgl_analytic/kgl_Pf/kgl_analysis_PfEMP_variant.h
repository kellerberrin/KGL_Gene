//
// Created by kellerberrin on 9/04/23.
//

#ifndef KGL_ANALYSIS_PFEMP_VARIANT_H
#define KGL_ANALYSIS_PFEMP_VARIANT_H


#include "kgl_genome_collection.h"
#include "kgl_variant_db_population.h"



namespace kellerberrin::genome {   //  organization::project level namespace



// A results class that holds the results of gene variant analysis for all genomes.
class GenomeGeneVariantAnalysis {

  using VariantGeneMap = std::map<std::string, std::pair<std::shared_ptr<const GeneFeature>, std::shared_ptr<ContigDB>>>;

public:

  GenomeGeneVariantAnalysis() = default;
  ~GenomeGeneVariantAnalysis() = default;

  void getGeneVariants(const std::shared_ptr<const PopulationDB> &population_ptr);
  void writeGeneResults(const std::string &file_name);
  void createGeneMap(const GeneVector& gene_vector);

private:

  VariantGeneMap variant_gene_map_;   // Population variant of selected genes.

  constexpr static const char CSV_DELIMITER_ = ',';

};



}; // Namespace.


#endif // KGL_ANALYSIS_PFEMP_VARIANT_H
