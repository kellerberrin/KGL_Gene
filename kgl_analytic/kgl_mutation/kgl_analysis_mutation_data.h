//
// Created by kellerberrin on 30/6/21.
//

#ifndef KGL_ANALYSIS_MUTATION_DATA_H
#define KGL_ANALYSIS_MUTATION_DATA_H

#include <vector>
#include <string>



namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convenience object. Somewhere to put lists of genes so that these don't have to be retrieved from a file.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class MutationAnalysisData {

public:

  MutationAnalysisData() = delete;
  ~MutationAnalysisData() = delete;

  // Prospective malaria active genes generated using GO Ontology similarity.
  [[nodiscard]] static const std::vector<std::string> &OntologyGeneList() { return ontology_derived_gene_list_; }

  // Malaria active genes from the OMIM entry #611162 available at https://www.omim.org/entry/611162
  [[nodiscard]] static const std::vector<std::string> OMIMGeneSymbol();
  [[nodiscard]] static const std::vector<std::string> OMIMGeneEnsembl();


  // Malaria active genes harvested from the Uniprot website.
  [[nodiscard]] static const std::vector<std::string> UniprotGeneSymbol();
  [[nodiscard]] static const std::vector<std::string> UniprotGeneEnsembl();

  // Used in adhoc gene polymorphism analysis.
  [[nodiscard]] static const std::vector<std::string> adHocGenes() { return adhoc_ensembl_symbol_; }

private:

  // From the OMIM entry #611162 available at https://www.omim.org/entry/611162
  static const std::vector<std::pair<std::string, std::string>> omim_ensembl_symbol_;

  // From the Uniprot website.
  static const std::vector<std::pair<std::string, std::string>> uniprot_ensembl_symbol_;

  // The list of genes to be analyzed variant by variant. Must be ensembl codes (for now).
  static const std::vector<std::string> ontology_derived_gene_list_;

  // Used in adhoc gene polymorphism analysis.
  static const std::vector<std::string> adhoc_ensembl_symbol_;

};


} //namespace

#endif //KGL_ANALYSIS_MUTATION_DATA_H
