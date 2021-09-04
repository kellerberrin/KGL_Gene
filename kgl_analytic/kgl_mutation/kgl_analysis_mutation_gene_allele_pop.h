//
// Created by kellerberrin on 31/8/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ALLELE_POP_H
#define KGL_ANALYSIS_MUTATION_GENE_ALLELE_POP_H


#include "kgl_uniprot_parser.h"
#include "kgl_entrez_parser.h"
#include "kgl_citation_parser.h"
#include "kgl_pubmed_resource.h"
#include "kgl_hsgenome_aux.h"
#include "kgl_variant_sort_analysis.h"
#include "kgl_variant_db_population.h"
#include "kgl_analysis_mutation_gene_ethnic.h"



namespace kellerberrin::genome {   //  organization::project level namespace


class GeneratePopulationAllele {

  using VariantCountMap = std::map<std::string, std::set<std::string>>;
  using ThreadReturnType = std::pair<std::string, std::set<std::string>>;

public:

  GeneratePopulationAllele() = default;
  ~GeneratePopulationAllele() = default;

  void initialize(const std::shared_ptr<const HsGenomeAux>& genome_aux_ptr,
                  const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                  const std::shared_ptr<const EntrezResource>& entrez_nomenclature_ptr,
                  const std::shared_ptr<const CitationResource>& allele_citation_ptr,
                  const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr);

  void processPopulation(const std::shared_ptr<const PopulationDB>& population_ptr);
  void writePopLiterature(const std::string& output_file) const;
  void addDiseaseAlleles(const DBCitationMap& disease_allele_map);

private:

  std::shared_ptr<const HsGenomeAux> genome_aux_ptr_;
  std::shared_ptr<const UniprotResource> uniprot_nomenclature_ptr_;
  std::shared_ptr<const EntrezResource> entrez_nomenclature_ptr_;
  std::shared_ptr<const CitationResource> allele_citation_ptr_;
  std::shared_ptr<const PubmedRequester> pubmed_requestor_ptr_;

  DBCitationMap disease_allele_map_;
  VariantCountMap variant_allele_map_;
  GeneEthnicitySex reference_ethnic_;

  const static constexpr char CONCAT_DELIMITER_{'&'};

  [[nodiscard]] std::set<std::string> getDiseaseCitations(const std::string& rs_code) const;
  [[nodiscard]] std::pair<std::string, std::string> generateGeneCodes(const std::vector<std::string>& ensembl_entrez_codes) const;
  [[nodiscard]] static ThreadReturnType getGenomePublications( std::shared_ptr<const GenomeDB> genome_ptr,
                                                               std::shared_ptr<const DBCitationMap> disease_cited_alleles_ptr);

};



} // namespace



#endif // KGL_ANALYSIS_MUTATION_GENE_ALLELE_POP_H
