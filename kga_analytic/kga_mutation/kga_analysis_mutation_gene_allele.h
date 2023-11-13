//
// Created by kellerberrin on 20/6/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
#define KGL_ANALYSIS_MUTATION_GENE_ALLELE_H


#include "kga_analysis_mutation_gene.h"
#include "kgl_variant_sort_analysis.h"
#include "kgl_pubmed_resource.h"



namespace kellerberrin::genome::analysis {   //  organization::project level namespace



class GenerateGeneAllele {

public:

  GenerateGeneAllele() = default;
  ~GenerateGeneAllele() = default;

  void initialize(const std::vector<std::string>& symbol_gene_list,
                  const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                  const std::shared_ptr<const EntrezResource>& entrez_nomenclature_ptr,
                  const std::shared_ptr<const CitationResource>& allele_citation_ptr,
                  const std::shared_ptr<const PubmedRequester>& pubmed_requestor_ptr);

  void writeOutput(const std::string& output_file, char delimiter) const;
  void addGeneCitedVariants(const std::shared_ptr<const SortedVariantAnalysis>& sorted_variants);
  void addDiseaseCitedVariants(const std::shared_ptr<const SortedVariantAnalysis>& sorted_variants);
  void addDiseaseAlleles(DBCitationMap disease_allele_map) { disease_allele_map_ = std::move(disease_allele_map); }
  // For each publication list the alleles.
  void writeLiteratureAlleleSummary(const std::string& output_file);
  // For each Allele print all the relevant publications.
  void writeAlleleLiteratureSummary(const std::string& output_file);

  void clear() { cited_allele_map_.clear(); }

private:


  const static constexpr char *AFR_SUPER_POP_{"AFR"};
  const static constexpr char *ALL_SUPER_POP_{"ALL"};

  // The key is the allele unique 'rs' code, the value pair is .first is a pointer to the variant, .second is a vector of all the gene codes for the allele.
  using RSCodeKey = std::string;
  using GeneCodes = std::vector<std::string>;
  using VariantPointer = std::shared_ptr<const Variant>;
  using AllelePair = std::pair<VariantPointer, GeneCodes>;
  using AlleleMap = std::map<RSCodeKey, AllelePair>;
  AlleleMap cited_allele_map_;

  std::shared_ptr<const UniprotResource> uniprot_nomenclature_ptr_;
  std::shared_ptr<const EntrezResource> entrez_nomenclature_ptr_;
  std::shared_ptr<const CitationResource> allele_citation_ptr_;
  std::shared_ptr<const PubmedRequester> pubmed_requestor_ptr_;
  std::map<std::string, std::string> ensembl_symbol_map_;
  DBCitationMap disease_allele_map_;

  static void writeHeader(std::ofstream& outfile, char delimiter);

  const static constexpr char CONCATENATE_VEP_FIELDS_{'&'};
  inline const static std::vector<std::string> VEP_FIELD_LIST_{ "Amino_acids",
                                                                "APPRIS",
                                                                "BIOTYPE",
                                                                "CANONICAL",
                                                                "CCDS",
                                                                "CDS_position",
                                                                "Codons",
                                                                "Consequence",
                                                                "Feature",
                                                                "Feature_type",
                                                                "GENE_PHENO",
                                                                "IMPACT",
                                                                "LoF",
                                                                "MOTIF_NAME",
                                                                "MOTIF_POS",
                                                                "MOTIF_SCORE_CHANGE",
                                                                "Protein_position",
                                                                "TSL"};

  [[nodiscard]] static std::map<std::string, std::string> retrieveVepFields( const std::shared_ptr<const Variant>& variant_ptr,
                                                                             const std::vector<std::string>& field_list);
  [[nodiscard]] std::set<std::string> getCitations(const std::string& rs_code) const;
  [[nodiscard]] bool citationsExist(const std::string& rs_code) const;
  [[nodiscard]] std::set<std::string> getDiseaseCitations(const std::string& rs_code) const;
  [[nodiscard]] std::pair<std::string, std::string> generateGeneCodes(const std::vector<std::string>& ensembl_entrez_codes) const;


};



} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
