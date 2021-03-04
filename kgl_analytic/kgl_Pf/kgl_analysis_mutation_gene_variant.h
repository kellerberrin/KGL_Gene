//
// Created by kellerberrin on 4/3/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_VARIANT_H
#define KGL_ANALYSIS_MUTATION_GENE_VARIANT_H


#include "kgl_variant_db_population.h"



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


  size_t male_phase_lof{0};          // Loss of gene function in the (B) chromosome.
  size_t female_phase_lof{0};        // Los of function in the female_ (A) chromosome.
  size_t hom_lof{0};          // Loss of gene function in both chromosomes.
  size_t male_lof{0};
  size_t female_lof{0};
  size_t male_high_effect{0};
  size_t female_high_effect{0};
  size_t hom_high_effect{0};          // Loss of gene function in both chromosomes.
  size_t male_moderate_effect{0};
  size_t female_moderate_effect{0};
  size_t hom_moderate_effect{0};          // Loss of gene function in both chromosomes.
  size_t male_modifier_effect{0};
  size_t female_modifier_effect{0};
  size_t hom_modifier_effect{0};          // Loss of gene function in both chromosomes.

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

  void writeVariantHeader(std::ostream &out_file, char output_delimiter) const;

  bool writeVariantOutput(std::ostream &out_file, char output_delimiter) const;

  void ProcessVariantStats(const std::shared_ptr<const ContigDB> &span_variant_ptr,
                           const std::shared_ptr<const PopulationDB> &unphased_population_ptr);

private:

  size_t unique_variants_{0};
  size_t span_variant_count_{0};
  size_t variant_count_{0};
  size_t male_phase_{0};  // Variants from the male_ phased (B) chromosome.
  size_t female_phase_{0};  // Variants from the female_ phased (A) chromosome.
  size_t male_lof_{0};          // Loss of gene function in the (B) chromosome.
  size_t female_lof_{0};        // Los of function in the female_ (A) chromosome.
  size_t hom_lof_{0};          // Loss of gene function in both chromosomes.
  size_t male_high_effect_{0};          // High Variant Impact in the (B) chromosome.
  size_t female_high_effect_{0};        // High Variant Impact in the (A) chromosome.
  size_t hom_high_effect_{0};          // High Impact in both in both chromosomes.
  size_t genome_count_{0};   // Total number of genomes.
  size_t genome_variant_{0};  // Number of genomes that contain variants for this gene.
  size_t homozygous_{0};
  size_t heterozygous_{0};
  double indel_{0.0};
  double transition_{0.0};
  double transversion_{0.0};

  // Vep fields.
  constexpr static const char *LOF_VEP_FIELD_ = "LoF";
  constexpr static const char *LOF_HC_VALUE_ = "HC";
  constexpr static const char *IMPACT_VEP_FIELD_ = "IMPACT";
  constexpr static const char *IMPACT_MODERATE_VALUE_ = "MODERATE";
  constexpr static const char *IMPACT_HIGH_VALUE_ = "HIGH";

  VepInfo geneSpanVep(const std::shared_ptr<const ContigDB> &span_contig,
                      const std::shared_ptr<const PopulationDB> &unphased_population_ptr);

  size_t VepCount(const std::shared_ptr<const ContigDB> &vep_contig,
                  const std::string &vep_field_ident,
                  const std::string &vep_field_value);

};


} // namespace


#endif //KGL_ANALYSIS_MUTATION_GENE_VARIANT_H
