//
// Created by kellerberrin on 20/6/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
#define KGL_ANALYSIS_MUTATION_GENE_ALLELE_H


#include "kgl_analysis_mutation_gene.h"



namespace kellerberrin::genome {   //  organization::project level namespace



class GenerateGeneAllele {

public:

  GenerateGeneAllele(const std::vector<std::string>& ensembl_gene_list) : sorted_allele_map_(std::make_shared<EnsemblIndexMap>()),
                                                                          ensembl_gene_list_(ensembl_gene_list) {}
  ~GenerateGeneAllele() = default;

  void updateAlleleMap(std::shared_ptr<const PopulationDB> unphased_population_ptr);
  void filterAlleleMap(const double AFR_frequency, const double upper_tail, const double lower_tail);

  void writeOutput(const std::string& output_file, const char delimiter) const;


private:


  const static constexpr char *AFR_SUPER_POP_{"AFR"};
  const static constexpr char *ALL_SUPER_POP_{"ALL"};

  std::shared_ptr<EnsemblIndexMap> sorted_allele_map_;
  const std::vector<std::string> ensembl_gene_list_;

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
