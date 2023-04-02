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
// High level object extracts resource specific properties.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ResourceProperties {

public:

  ResourceProperties(std::shared_ptr<const PropertyTree> property_tree_ptr, std::string work_directory)
    : property_tree_ptr_(std::move(property_tree_ptr)), work_directory_(std::move(work_directory)) {}
  ~ResourceProperties() = default;

  [[nodiscard]]  ResourceDefinitions getRuntimeResources() const;

  // Resource Database categories.
  constexpr static const char RESOURCE_LIST_[] = "resourceList";
  // Go Ontology Database; this is a Resource Database.
  constexpr static const char ONTOLOGY_DATABASE_[] = "ontologyDatabase";
  constexpr static const char ONTOLOGY_IDENT_[] = "ontologyIdent";
  constexpr static const char GAF_ANNOTATION_FILE_[] = "gafFile";
  constexpr static const char GO_ONTOLOGY_FILE_[] = "goFile";
  // Genome Database categories; this a Resource Database.
  constexpr static const char GENOME_DATABASE_[] = "genomeDatabase";
  constexpr static const char GENOME_IDENT_[] = "genomeIdent";
  constexpr static const char FASTA_FILE_[] = "fastaFile";
  constexpr static const char GFF_FILE_[] = "gffFile";
  constexpr static const char TRANSLATION_TABLE_[] = "translationTable";
  // Gene Nomenclature files; this a Resource Database.
  constexpr static const char GENE_ID_DATABASE_[] = "geneNomenclature";
  constexpr static const char GENE_ID_IDENT_[] = "nomenclatureIdent";
  // Gene Nomenclature identifiers.
  static const constexpr char* NOMENCLATURE_UNIPROTID{"UniprotID"};   // The uniprot nomenclature file
  static const constexpr char* NOMENCLATURE_ENSEMBL{"EnsemblHGNC"};   // The ensembl nomenclature file
  constexpr static const char GENE_ID_FILE_[] = "nomenclatureFile";
  // genome genealogy files; this a Resource Database.
  constexpr static const char GENEALOGY_ID_DATABASE_[] = "genomeGenealogy";
  constexpr static const char GENEALOGY_ID_IDENT_[] = "genealogyIdent";
  constexpr static const char GENEALOGY_ID_FILE_[] = "genealogyFile";
  // genome aux info files; this a Resource Database.
  constexpr static const char AUX_ID_DATABASE_[] = "genomeAux";
  constexpr static const char AUX_ID_IDENT_[] = "auxIdent";
  constexpr static const char AUX_ID_FILE_[] = "auxFile";
  // Allele Citations; this is a resource database.
  constexpr static const char CITATION_DATABASE_[] = "alleleCitation";
  constexpr static const char CITATION_IDENT_[] = "citationIdent";
  constexpr static const char CITATION_FILE_[] = "citationFile";
  // Human Entrez Gene Infomation; this is a resource database.
  constexpr static const char ENTREZ_DATABASE_[] = "entrezGene";
  constexpr static const char ENTREZ_IDENT_[] = "entrezIdent";
  constexpr static const char ENTREZ_FILE_[] = "entrezFile";
  // Pf7 sample Infomation; this is a resource database.
  constexpr static const char PF7_SAMPLE_DATABASE_[] = "Pf7Sample";
  constexpr static const char PF7_SAMPLE_IDENT_[] = "Pf7SampleIdent";
  constexpr static const char PF7_SAMPLE_FILE_[] = "Pf7SampleFile";
  // Literature PMID linked to BioConcepts (Gene, Disease, etc); this is a resource database.
  constexpr static const char PMID_BIO_DATABASE_[] = "bioPMID";
  constexpr static const char PMID_BIO_IDENT_[] = "bioPMIDIdent";
  constexpr static const char PMID_BIO_FILE_[] = "bioPMIDFile";
  // The pubmed rest literature API; this is a resource api.
  constexpr static const char PUBMED_LIT_API_[] = "pubmedApi";
  constexpr static const char PUBMED_LIT_IDENT_[] = "pubmedApiIdent";
  constexpr static const char PUBMED_PUBLICATION_CACHE_[] = "pubmedPublicationCache";
  constexpr static const char PUBMED_CITATION_CACHE_[] = "pubmedCitationCache";

private:

  std::shared_ptr<const PropertyTree> property_tree_ptr_;   // The aggregated and parsed XML property tree.
  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'

  [[nodiscard]] const std::string& workDirectory() const { return work_directory_; }

  constexpr static const char DOT_[] = ".";
  constexpr static const char RUNTIME_ROOT_[] = "runTime";

  [[nodiscard]] std::optional<ResourceParameters> genomeResource(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> ontologyDatabase(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> geneIDDatabase(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> genealogyIDDatabase(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> citationDatabase(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> entrezDatabase(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> Pf7SampleDatabase(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> PubmedLitAPI(const PropertyTree& sub_tree) const;
  [[nodiscard]] std::optional<ResourceParameters> auxIDDatabase(const PropertyTree& sub_tree) const;

};


}   // end namespace


#endif //KGL_PROPERTIES_RESOURCE_H
