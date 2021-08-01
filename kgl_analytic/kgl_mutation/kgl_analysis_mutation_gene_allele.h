//
// Created by kellerberrin on 20/6/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
#define KGL_ANALYSIS_MUTATION_GENE_ALLELE_H


#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_sort_analysis.h"



namespace kellerberrin::genome {   //  organization::project level namespace



class GenerateGeneAllele {

public:

  GenerateGeneAllele() = default;
  ~GenerateGeneAllele() = default;

  void initialize(const std::vector<std::string>& symbol_gene_list,
                  const std::shared_ptr<const UniprotResource>& uniprot_nomenclature_ptr,
                  const std::shared_ptr<const CitationResource>& allele_citation_ptr);

  void filterAlleleMap(const double ALL_frequency, const double AFR_frequency);
  void writeOutput(const std::string& output_file, const char delimiter) const;
  void addSortedVariants(const std::shared_ptr<const SortedVariantAnalysis>& sorted_variants);

private:


  const static constexpr char *AFR_SUPER_POP_{"AFR"};
  const static constexpr char *ALL_SUPER_POP_{"ALL"};

  EnsemblIndexMap sorted_allele_map_;
  std::shared_ptr<const CitationResource> allele_citation_ptr_;
  std::map<std::string, std::string> ensembl_symbol_map_;

  static void writeHeader(std::ofstream& outfile, const char delimiter);

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

  static std::map<std::string, std::string> retrieveVepFields( const std::shared_ptr<const Variant>& variant_ptr,
                                                               const std::vector<std::string>& field_list);

};



} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
