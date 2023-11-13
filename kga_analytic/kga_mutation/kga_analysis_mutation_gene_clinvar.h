//
// Created by kellerberrin on 31/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H
#define KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H

#include "kgl_variant_db_population.h"
#include "kga_analysis_mutation_gene_ethnic.h"



namespace kellerberrin::genome::analysis {   //  organization::project level namespace


struct ClinvarInfo {

  std::string clndn;         // Clinical description.
  std::shared_ptr<const Variant> variant_ptr;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneClinvar {

public:

  GeneClinvar()   {

    clinvar_contig_ = std::make_shared<ContigDB>("");
    clinvar_ethnic_.setDisplay("ETH_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));

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
                      const ContigId_t& contig_id,
                      const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                      const std::shared_ptr<const ContigDB>& gene_variants,
                      const std::shared_ptr<const HsGenomeAux>& genome_aux_data);


  // Superpopulation, population and sex breakdown.
  [[nodiscard]] GeneEthnicitySex& updateEthnicity() { return clinvar_ethnic_; }


private:

  inline constinit const static char* CONCAT_TOKEN_ = "&";
  // Text description of problems with this gene.
  std::set<std::string> clinvar_desc_;
  // Superpopulation and sex breakdown.
  GeneEthnicitySex clinvar_ethnic_;
  size_t hom_genome_{0};  // Homozygous
  size_t genome_count_{0};
  std::shared_ptr<const ContigDB> clinvar_contig_;

  // Clinvar fields.
  constexpr static const char* CLINVAR_CLNDN_FIELD = "CLNDN";
  constexpr static const char* CLINVAR_CLNSIG_FIELD = "CLNSIG";
  constexpr static const char* CLINVAR_PATH_SIGNIF = "PATHOGENIC";

  // Text description of problems with this gene.
  [[nodiscard]] const std::set<std::string>& getClinvarDesc() const { return  clinvar_desc_; }

  // Superpopulation, population and sex breakdown.
  [[nodiscard]] const GeneEthnicitySex& getEthnicity() const { return clinvar_ethnic_; }

  void processClinvar(const GenomeId_t& genome_id,
                      const std::shared_ptr<const ContigDB>& gene_variants,
                      const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  static std::vector<ClinvarInfo> clinvarInfo(const std::shared_ptr<const ContigDB>& clinvar_contig_ptr);

  static std::shared_ptr<const ContigDB> getClinvarContig(const ContigId_t& contig_id,
                                                          const std::shared_ptr<const PopulationDB>& clinvar_population_ptr);
  static std::shared_ptr<const ContigDB> FilterPathogenic(std::shared_ptr<const ContigDB> clinvar_contig);



};


} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H
