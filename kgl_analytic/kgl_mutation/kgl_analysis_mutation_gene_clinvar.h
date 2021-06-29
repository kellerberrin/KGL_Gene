//
// Created by kellerberrin on 31/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H
#define KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H

#include "kgl_variant_db_population.h"
#include "kgl_analysis_mutation_gene_ethnic.h"



namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantPhaseStats {

public:

  VariantPhaseStats() = default;
  VariantPhaseStats(const VariantPhaseStats &) = default;
  ~VariantPhaseStats() = default;

  VariantPhaseStats &operator=(const VariantPhaseStats &) = default;

  bool phaseStatistics(const std::shared_ptr<const ContigDB>& contig_ptr);
  static void writeHeader(std::ostream& out_file, char output_delimiter);

  [[nodiscard]] size_t phaseMale() const { return phase_male_; }  // Allele(s) is present on the male_ phase
  [[nodiscard]] size_t phaseFemale() const { return phase_female_; }  // Allele(s) is present on the female_ phase
  [[nodiscard]] size_t phaseHomozygous() const { return phase_hom_; }  // Allele(s) is present on both phases. Male and female_.
  [[nodiscard]] size_t phaseEither() const { return phase_either_; }  // Allele(s) is present on either phase. Male or female_.

  [[nodiscard]] size_t& updatePhaseMale() { return phase_male_; }  // Allele(s) is present on the male_ phase
  [[nodiscard]] size_t& updatePhaseFemale() { return phase_female_; }  // Allele(s) is present on the female_ phase
  [[nodiscard]] size_t& updatePhaseHomozygous() { return phase_hom_; }  // Allele(s) is present on both phases. Male and female_.
  [[nodiscard]] size_t& updatePhaseEither() { return phase_either_; }  // Allele(s) is present on either phase. Male or female_.

private:

  // Variant phase
  size_t phase_male_{0};   // Allele(s) is present on the male_ phase
  size_t phase_female_{0};  // Allele(s) is present on the female_ phase
  size_t phase_hom_{0};  // Allele(s) is present on both phases. Male and female_.
  size_t phase_either_{0};  // Allele(s) is present on either phase. Male or female_.

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneClinvar {

public:

  GeneClinvar() {

    results_.setDisplay("ETH_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));

  }

  GeneClinvar(const GeneClinvar &) = default;

  ~GeneClinvar() = default;

  GeneClinvar &operator=(const GeneClinvar &) = default;


  void writeOutput( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                    std::ostream& out_file,
                    char output_delimiter) const;

  void writeHeader( const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                    std::ostream& out_file,
                    char output_delimiter) const;

  void processClinvar(const GenomeId_t& genome_id,
                      const std::shared_ptr<const ContigDB>& gene_variants,
                      const std::shared_ptr<const ContigDB>& clinvar_contig,
                      const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  // Superpopulation, population and sex breakdown.
  [[nodiscard]] GeneEthnicitySex& updateEthnicity() { return results_; }

private:

  inline constinit const static char* CONCAT_TOKEN_ = "&";
  // Variant phase
  VariantPhaseStats phase_;
  // Text description of problems with this gene.
  std::set<std::string> clinvar_desc_;
  // Superpopulation and sex breakdown.
  GeneEthnicitySex results_;

  // Variant phase information.
  [[nodiscard]] const VariantPhaseStats& getPhase() const { return phase_; }

  // Text description of problems with this gene.
  [[nodiscard]] const std::set<std::string>& getClinvarDesc() const { return  clinvar_desc_; }

  // Superpopulation, population and sex breakdown.
  [[nodiscard]] const GeneEthnicitySex& getEthnicity() const { return results_; }

  // Variant phase information.
  [[nodiscard]] VariantPhaseStats& updatePhase() { return phase_; }

  // Text description of problems with this gene.
  [[nodiscard]] std::set<std::string>& updateClinvarDesc() { return  clinvar_desc_; }



};


} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H
