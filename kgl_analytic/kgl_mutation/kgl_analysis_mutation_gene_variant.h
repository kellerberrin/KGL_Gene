//
// Created by kellerberrin on 4/3/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_VARIANT_H
#define KGL_ANALYSIS_MUTATION_GENE_VARIANT_H

#include "kgl_citation_parser.h"
#include "kgl_variant_db_population.h"
#include "kgl_analysis_mutation_gene_ethnic.h"



namespace kellerberrin::genome {   //  organization::project level namespace




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VepInfo {

public:

  VepInfo() = default;

  ~VepInfo() = default;

  VepInfo(const VepInfo &) = default;

  VepInfo &operator=(const VepInfo &) = default;


  size_t all_lof{0};
  size_t hom_lof{0};          // Loss of gene function in both chromosomes.
  size_t all_high_effect{0};
  size_t hom_high_effect{0};          // Loss of gene function in both chromosomes.
  size_t all_moderate_effect{0};
  size_t hom_moderate_effect{0};          // Loss of gene function in both chromosomes.

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneVariants {

public:

  GeneVariants() = default;
  ~GeneVariants() = default;

  GeneVariants(const GeneVariants &) = default;

  GeneVariants &operator=(const GeneVariants &) = default;

  // Must be called before updating object.
  void initializeEthnic(const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  void writeVariantHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                           std::ostream &out_file,
                           char output_delimiter) const;

  bool writeVariantOutput( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                           std::ostream &out_file,
                           char output_delimiter) const;

  void processVariantStats(const GenomeId_t& genome,
                           const std::shared_ptr<const ContigDB> &span_variant_ptr,
                           const std::shared_ptr<const PopulationDB> &unphased_population_ptr,
                           const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                           const std::shared_ptr<const CitationResource>& allele_citation_ptr);

  [[nodiscard]] bool processSummaryStatistics( const std::shared_ptr<const PopulationDB> &population_ptr,
                                               const GeneEthnicitySex& ethnic_statistics,
                                               const std::string& gene);

  void initializeSummaryStatistics( const GeneEthnicitySex& ethnic_statistics);


private:

  std::set<std::string> unique_variants_;
  size_t span_variant_count_{0};
  size_t variant_count_{0};
  size_t all_lof_{0};            // All lof variants;
  size_t hom_lof_{0};
  GeneEthnicitySex ethnic_lof_;
  size_t all_high_effect_{0};          // High Variant Impact
  size_t hom_high_effect_{0};
  GeneEthnicitySex ethnic_high_;
  size_t all_moderate_effect_{0};          // Moderate impact variant
  size_t hom_moderate_effect_{0};
  GeneEthnicitySex ethnic_moderate_;
  size_t citation_count_{0};
  std::set<std::string> unique_citations_;
  GeneEthnicitySex ethnic_citation_;
  size_t genome_count_{0};   // Total number of genomes.
  size_t genome_variant_{0};  // Number of genomes that contain variants for this gene.

  // Hypergeometric summary statistics
  std::map<std::string, double> upper_tail_;
  std::map<std::string, double> lower_tail_;

  // Vep fields.
  constexpr static const char *LOF_VEP_FIELD_ = "LoF";
  constexpr static const char *LOF_HC_VALUE_ = "HC";
  constexpr static const char *IMPACT_VEP_FIELD_ = "IMPACT";
  constexpr static const char *IMPACT_MODERATE_VALUE_ = "MODERATE";
  constexpr static const char *IMPACT_HIGH_VALUE_ = "HIGH";

  VepInfo geneSpanVep(const std::shared_ptr<const ContigDB> &span_contig,
                      const std::shared_ptr<const PopulationDB> &unphased_population_ptr);

  size_t vepCount(const std::shared_ptr<const ContigDB> &vep_contig,
                  const std::string &vep_field_ident,
                  const std::string &vep_field_value);

};


} // namespace


#endif //KGL_ANALYSIS_MUTATION_GENE_VARIANT_H
