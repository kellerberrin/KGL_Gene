//
// Created by kellerberrin on 20/6/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
#define KGL_ANALYSIS_MUTATION_GENE_ALLELE_H


#include "kgl_analysis_mutation_gene.h"



namespace kellerberrin::genome {   //  organization::project level namespace



class GenerateGeneAllele {

public:

  GenerateGeneAllele() : sorted_allele_map_(std::make_shared<EnsemblIndexMap>()) {}
  ~GenerateGeneAllele() = default;

  void updateAlleleMap(std::shared_ptr<const PopulationDB> unphased_population_ptr);
  void filterAlleleMap(const double AFR_frequency, const double upper_tail, const double lower_tail);

  void writeOutput(const std::string& output_file, const char delimiter) const;


private:


  const static constexpr char *AFR_SUPER_POP_{"AFR"};
  const static constexpr char *ALL_SUPER_POP_{"ALL"};

  const std::vector<std::string> ensembl_gene_list_{

      "ENSG00000284690", "ENSG00000282992", "ENSG00000275019", "ENSG00000262576", "ENSG00000256797",
      "ENSG00000254521", "ENSG00000213402", "ENSG00000204345", "ENSG00000198178", "ENSG00000196371",
      "ENSG00000189184", "ENSG00000188211", "ENSG00000186407", "ENSG00000185475", "ENSG00000185187",
      "ENSG00000183840", "ENSG00000183019", "ENSG00000180549", "ENSG00000179213", "ENSG00000172794",
      "ENSG00000172322", "ENSG00000171840", "ENSG00000170956", "ENSG00000170425", "ENSG00000169704",
      "ENSG00000167850", "ENSG00000167123", "ENSG00000166589", "ENSG00000165682", "ENSG00000164713",
      "ENSG00000163600", "ENSG00000163485", "ENSG00000162897", "ENSG00000161649", "ENSG00000159674"

  };

  std::shared_ptr<EnsemblIndexMap> sorted_allele_map_;

  static void writeHeader(std::ofstream& outfile, const char delimiter);


};



} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_ALLELE_H
