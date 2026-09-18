//
// Created by kellerberrin on 2/04/23.
//

#ifndef KGL_PROPERTIES_RESOURCE_H
#define KGL_PROPERTIES_RESOURCE_H


#include "kgl_runtime_resource.h"

#include <string>


namespace kellerberrin::genome {   //  organization::project level namespace


/// Compatibility shim exposing the legacy resource identifier constants.
/// The canonical definitions now live in ResourceType (kgl_resource_type.h).
/// This header exists so that code outside kgl_app (kgl_genomics, kol_ontology,
/// kga_analytic) that references ResourceProperties::*_RESOURCE_ID_ continues to compile.
class ResourceProperties {

public:

  /// XML element name for the resource list.
  constexpr static const char RESOURCE_LIST_[] = "resourceList";

  /// Genome database resource type identifier.
  constexpr static const char GENOME_RESOURCE_ID_[] = "genomeDatabase";
  /// Ontology database resource type identifier.
  constexpr static const char ONTOLOGY_RESOURCE_ID_[] = "ontologyDatabase";
  /// Gene nomenclature resource type identifier.
  constexpr static const char GENE_NOMENCLATURE_RESOURCE_ID_[] = "geneNomenclature";
  /// Genome genealogy resource type identifier.
  constexpr static const char GENEALOGY_RESOURCE_ID_[] = "genomeGenealogy";
  /// Genome auxiliary information resource type identifier.
  constexpr static const char GENOMEAUX_RESOURCE_ID_[] = "genomeAux";
  /// Allele citation resource type identifier.
  constexpr static const char CITATION_RESOURCE_ID_[] = "alleleCitation";
  /// Entrez gene resource type identifier.
  constexpr static const char ENTREZ_RESOURCE_ID_[] = "entrezGene";
  /// PubMed API resource type identifier.
  constexpr static const char PUBMED_API_RESOURCE_ID_[] = "pubmedApi";
  /// Pf7 sample data resource type identifier.
  constexpr static const char PF7SAMPLE_RESOURCE_ID_[] = "Pf7Sample";
  /// Pf7 FWS resource type identifier.
  constexpr static const char PF7FWS_RESOURCE_ID_[] = "Pf7Fws";
  /// Pf7 genetic distance resource type identifier.
  constexpr static const char PF7DISTANCE_RESOURCE_ID_[] = "Pf7Distance";
  /// Pf3K COI resource type identifier.
  constexpr static const char PF3K_COI_RESOURCE_ID_[] = "Pf3KCOI";

  /// Gene nomenclature identifier for Uniprot ID.
  static const constexpr char* NOMENCLATURE_UNIPROTID{"UniprotID"};
  /// Gene nomenclature identifier for Ensembl HGNC.
  static const constexpr char* NOMENCLATURE_ENSEMBL{"EnsemblHGNC"};

};


}   // end namespace


#endif //KGL_PROPERTIES_RESOURCE_H