//
// Created by kellerberrin on 4/3/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_VARIANT_H
#define KGL_ANALYSIS_MUTATION_GENE_VARIANT_H


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

  GeneVariants() {

    ethnic_lof_.setDisplay("LOF_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));
    ethnic_high_.setDisplay("HIGH_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));
    ethnic_moderate_.setDisplay("MOD_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));

  }

  ~GeneVariants() = default;

  GeneVariants(const GeneVariants &) = default;

  GeneVariants &operator=(const GeneVariants &) = default;

  void writeVariantHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                           std::ostream &out_file,
                           char output_delimiter) const;

  bool writeVariantOutput( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                           std::ostream &out_file,
                           char output_delimiter) const;

  void processVariantStats(const GenomeId_t& genome,
                           const std::shared_ptr<const ContigDB> &span_variant_ptr,
                           const std::shared_ptr<const PopulationDB> &unphased_population_ptr,
                           const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  [[nodiscard]] bool processSummaryStatistics( const std::shared_ptr<const PopulationDB> &population_ptr,
                                               const GeneEthnicitySex& ethnic_statistics,
                                               const std::string& gene);

  void initializeSummaryStatistics( const GeneEthnicitySex& ethnic_statistics);

  [[nodiscard]] GeneEthnicitySex& updateLofEthnicity() { return ethnic_lof_; }
  [[nodiscard]] GeneEthnicitySex& updateHighEthnicity() { return ethnic_high_; }
  [[nodiscard]] GeneEthnicitySex& updateModerateEthnicity() { return ethnic_moderate_; }

private:

  size_t unique_variants_{0};
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
