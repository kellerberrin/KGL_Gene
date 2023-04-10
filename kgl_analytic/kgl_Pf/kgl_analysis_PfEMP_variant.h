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

  GenomeGeneVariantAnalysis() {

    gene_population_ptr = std::make_unique<PopulationDB>("GenePopulation", DataSourceEnum::Falciparum);

  }
  ~GenomeGeneVariantAnalysis() = default;

  void writeGeneResults(const std::string &file_name);
  void setGeneVector(const GeneVector& gene_vector);
  void getGeneVariants(const std::shared_ptr<const PopulationDB> &population_ptr);

private:

  GeneVector gene_vector_;   // Gene vector of interest.
  std::unique_ptr<PopulationDB> gene_population_ptr; // Each contig per genome is a gene in the gene vector

  constexpr static const char CSV_DELIMITER_ = ',';



};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneGenomeAnalysis {

  using VariantGenomeCount = std::map<std::string, std::pair<std::shared_ptr<const Variant>, std::vector<GenomeId_t>>>;

public:

  explicit GeneGenomeAnalysis(const std::shared_ptr<const ContigDB>& gene_unique_variants);
  ~GeneGenomeAnalysis() = default;

private:

  VariantGenomeCount gene_genome_analysis_;

  constexpr static const char* NULL_VARIANT_ = "NullVariant";

};




}; // Namespace.


#endif // KGL_ANALYSIS_PFEMP_VARIANT_H
