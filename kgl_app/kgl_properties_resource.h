//
// Created by kellerberrin on 2/04/23.
//

#ifndef KGL_PROPERTIES_RESOURCE_H
#define KGL_PROPERTIES_RESOURCE_H


#include "kgl_runtime_resource.h"
#include "kel_property_tree.h"

#include <memory>
#include <string>
#include <map>
#include <set>


namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// High level object extracts resource specific properties from the XML tree.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Steps for adding a resource.
//
// 1. Define the resource in a resource XML file.
// 2. Define all the tags used in the XML definition below as text constants.
// 3. Go to the kgl_properties_resource.cpp and kgl_properties_resource_item.cpp and define a simple XML parser.
// 4. Go to kgl_package_resource.cpp function ExecutePackage::loadRuntimeResource() and add code to create the resource on demand.
// 5. Add the resource type and identifier in the package XML.
// 6. Retrieve the resource at runtime in the analysis object.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ResourceProperties {

public:

  ResourceProperties(std::shared_ptr<const PropertyTree> property_tree_ptr, std::string work_directory)
    : property_tree_ptr_(std::move(property_tree_ptr)), work_directory_(std::move(work_directory)) {}
  ~ResourceProperties() = default;

  [[nodiscard]]  ResourceDefinitions getRuntimeResources() const;

  // Resource Database categories.
  constexpr static const char RESOURCE_LIST_[] = "resourceList";
  // Go Ontology Database; this is a Resource Database.
  constexpr static const char ONTOLOGY_RESOURCE_ID_[] = "ontologyResourceID";
  constexpr static const char ONTOLOGY_IDENT_[] = "ontologyIdent";
  constexpr static const char GAF_ANNOTATION_FILE_[] = "gafFile";
  constexpr static const char GO_ONTOLOGY_FILE_[] = "goFile";
  // Genome Database categories; this a Resource Database.
  constexpr static const char GENOME_RESOURCE_ID_[] = "genomeDatabase";
  constexpr static const char GENOME_IDENT_[] = "genomeIdent";
  constexpr static const char FASTA_FILE_[] = "fastaFile";
  constexpr static const char GFF_FILE_[] = "gffFile";
  constexpr static const char TRANSLATION_TABLE_[] = "translationTable";
  // Gene Nomenclature files; this a Resource Database.
  constexpr static const char GENE_NOMENCLATURE_RESOURCE_ID_[] = "geneNomenclature";
  constexpr static const char GENE_NOMENCLATURE_IDENT_[] = "nomenclatureIdent";
  // Gene Nomenclature identifiers.
  static const constexpr char* NOMENCLATURE_UNIPROTID{"UniprotID"};   // The uniprot nomenclature file
  static const constexpr char* NOMENCLATURE_ENSEMBL{"EnsemblHGNC"};   // The ensembl nomenclature file
  constexpr static const char GENE_NOMENCLATURE_FILE_[] = "nomenclatureFile";
  // genome genealogy files; this a Resource Database.
  constexpr static const char GENEALOGY_RESOURCE_ID_[] = "genomeGenealogy";
  constexpr static const char GENEALOGY_IDENT_[] = "genealogyIdent";
  constexpr static const char GENEALOGY_FILE_[] = "genealogyFile";
  // genome aux info files; this a Resource Database.
  constexpr static const char GENOMEAUX_RESOURCE_ID_[] = "genomeAux";
  constexpr static const char GENOMEAUX_IDENT_[] = "auxIdent";
  constexpr static const char GENOMEAUX_FILE_[] = "auxFile";
  // Allele Citations; this is a resource database.
  constexpr static const char CITATION_RESOURCE_ID_[] = "alleleCitation";
  constexpr static const char CITATION_IDENT_[] = "citationIdent";
  constexpr static const char CITATION_FILE_[] = "citationFile";
  // Human Entrez Gene Infomation; this is a resource database.
  constexpr static const char ENTREZ_RESOURCE_ID_[] = "entrezGene";
  constexpr static const char ENTREZ_IDENT_[] = "entrezIdent";
  constexpr static const char ENTREZ_FILE_[] = "entrezFile";
  // Pf7 sample Infomation; this is a resource database.
  constexpr static const char PF7SAMPLE_RESOURCE_ID_[] = "Pf7Sample";
  constexpr static const char PF7SAMPLE_IDENT_[] = "Pf7SampleIdent";
  constexpr static const char PF7SAMPLE_FILE_[] = "Pf7SampleFile";
  // Pf7 within-host infection fixation index (FWS); this is a resource database.
  constexpr static const char PF7FWS_RESOURCE_ID_[] = "Pf7Fws";
  constexpr static const char PF7FWS_IDENT_[] = "Pf7FwsIdent";
  constexpr static const char PF7FWS_FILE_[] = "Pf7FwsFile";
  // Literature PMID linked to BioConcepts (Gene, Disease, etc); this is a resource database.
  constexpr static const char PMIDBIO_RESOURCE_ID_[] = "bioPMID";
  constexpr static const char PMIDBIO_IDENT_[] = "bioPMIDIdent";
  constexpr static const char PMIDBIO_FILE_[] = "bioPMIDFile";
  // The pubmed rest literature API; this is a resource api.
  constexpr static const char PUBMED_API_RESOURCE_ID_[] = "pubmedApi";
  constexpr static const char PUBMED_IDENT_[] = "pubmedApiIdent";
  constexpr static const char PUBMED_PUBLICATION_CACHE_[] = "pubmedPublicationCache";
  constexpr static const char PUBMED_CITATION_CACHE_[] = "pubmedCitationCache";

private:

  std::shared_ptr<const PropertyTree> property_tree_ptr_;   // The aggregated and parsed XML property tree.
  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'

  [[nodiscard]] const std::string& workDirectory() const { return work_directory_; }

  constexpr static const char DOT_[] = ".";
  constexpr static const char RUNTIME_ROOT_[] = "runTime";

  [[nodiscard]] std::optional<ResourceParameters> genomeResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> ontologyResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> geneIDResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> genealogyIDResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> citationResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> entrezResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> Pf7SampleResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> Pf7FwsResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> PubmedLitAPIResourceXML(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> auxIDResourceXML(const PropertyTree& sub_tree) const;

};


}   // end namespace


#endif //KGL_PROPERTIES_RESOURCE_H
