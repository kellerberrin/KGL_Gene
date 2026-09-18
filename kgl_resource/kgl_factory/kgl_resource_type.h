//
// Created by kellerberrin on 26/08/26.
//

#ifndef KGL_RESOURCE_TYPE_H
#define KGL_RESOURCE_TYPE_H

#include <string_view>

namespace kellerberrin::genome {   //  organization::project level namespace


/// Standalone resource type identifiers. These are used both as XML element names
/// under <resourceList> and as ResourceBase::resourceType() values.
struct ResourceType {

  /// Genome database resource identifier.
  static constexpr std::string_view GENOME{"genomeDatabase"};
  /// Ontology database resource identifier.
  static constexpr std::string_view ONTOLOGY{"ontologyDatabase"};
  /// Gene nomenclature resource identifier.
  static constexpr std::string_view GENE_NOMENCLATURE{"geneNomenclature"};
  /// Genome genealogy resource identifier.
  static constexpr std::string_view GENEALOGY{"genomeGenealogy"};
  /// Genome auxiliary information resource identifier.
  static constexpr std::string_view GENOME_AUX{"genomeAux"};
  /// Allele citation resource identifier.
  static constexpr std::string_view CITATION{"alleleCitation"};
  /// Entrez gene resource identifier.
  static constexpr std::string_view ENTREZ{"entrezGene"};
  /// PubMed API resource identifier.
  static constexpr std::string_view PUBMED_API{"pubmedApi"};
  /// Pf7 sample data resource identifier.
  static constexpr std::string_view PF7_SAMPLE{"Pf7Sample"};
  /// Pf7 FWS resource identifier.
  static constexpr std::string_view PF7_FWS{"Pf7Fws"};
  /// Pf7 genetic distance resource identifier.
  static constexpr std::string_view PF7_DISTANCE{"Pf7Distance"};
  /// Pf3K COI resource identifier.
  static constexpr std::string_view PF3K_COI{"Pf3KCOI"};

  /// Gene nomenclature sub-identifier for Uniprot ID selection at load time.
  static constexpr std::string_view NOMENCLATURE_UNIPROT{"UniprotID"};
  /// Gene nomenclature sub-identifier for Ensembl HGNC selection at load time.
  static constexpr std::string_view NOMENCLATURE_ENSEMBL{"EnsemblHGNC"};

};


}   // end namespace


#endif //KGL_RESOURCE_TYPE_H