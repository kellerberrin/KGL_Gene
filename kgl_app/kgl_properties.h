//
// Created by kellerberrin on 11/11/18.
//


#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H


#include "kgl_runtime.h"
#include "kel_property_tree.h"
#include "kgl_genome_types.h"

#include <memory>
#include <string>
#include <map>
#include <set>


namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// High level object extracts application specific properties.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeProperties {

public:

  RuntimeProperties() = default;
  ~RuntimeProperties() = default;

  [[nodiscard]] bool readProperties(const std::string& properties_file);

  void setWorkDirectory(const std::string& work_directory) { work_directory_ = work_directory; }

  [[nodiscard]] const std::string& workDirectory() const { return work_directory_; }

  [[nodiscard]] ActivePackageVector getActivePackages() const;

  [[nodiscard]] RuntimePackageMap getPackageMap() const;

  [[nodiscard]]  RuntimeAnalysisMap getAnalysisMap() const;

  [[nodiscard]]  ResourceDefinitions getRuntimeResourceDef() const;

  [[nodiscard]] RuntimeDataFileMap getDataFiles() const;

  [[nodiscard]] ContigAliasMap getContigAlias() const;

  [[nodiscard]] VariantEvidenceMap getEvidenceMap() const;

  [[nodiscard]] ActiveParameterList getParameterMap() const;

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
  // Ignore these.
  constexpr static const char HELP_[] = "help";
  constexpr static const char COMMENT_[] = "<xmlcomment>";

private:

  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'
  PropertyTree property_tree_;   // All the option XML files.

  // Node categories.
  constexpr static const char DOT_[] = ".";
  constexpr static const char RUNTIME_ROOT_[] = "runTime";
  constexpr static const char ACTIVE_[] = "active";
  constexpr static const char VALUE_[] = "value";

  // Active Package Runtime categories.
  constexpr static const char EXECUTE_LIST_[] = "executeList";
  // Package Runtime categories.
  constexpr static const char PACKAGE_LIST_[] = "packageList";
  constexpr static const char PACKAGE_[] = "package";
  constexpr static const char PACKAGE_IDENT_[] = "packageIdent";
  constexpr static const char PACKAGE_ANALYSIS_LIST_[] = "analysisList";
  constexpr static const char PACKAGE_RESOURCE_LIST_[] = "resourceList";
  constexpr static const char PACKAGE_ITERATION_[] = "iteration";
  constexpr static const char PACKAGE_ITERATION_LIST_[] = "iterationList";
  // Analysis Runtime categories.
  constexpr static const char ANALYSIS_LIST_[] = "analysisList";
  constexpr static const char ANALYSIS_[] = "analysis";
  constexpr static const char ANALYSIS_IDENT_[] = "analysisIdent";
  // Parameter Runtime categories.
  constexpr static const char PARAMETER_LIST_[] = "parameterList";
  constexpr static const char PARAMETER_BLOCK_[] = "parameterBlock";
  constexpr static const char PARAMETER_NAME_[] = "parameterName";
  constexpr static const char PARAMETER_VECTOR_[] = "parameterVector";
  constexpr static const char PARAMETER_[] = "parameter";
  constexpr static const char PARAMETER_IDENT_[] = "parameterIdent";
  constexpr static const char PARAMETER_VALUE_[] = "parameterValue";
  constexpr static const char PARAMETER_RUNTIME_[] = "parameterRuntime";
  // Data File Runtime categories.
  constexpr static const char DATA_FILE_LIST_[] = "dataFileList";
  constexpr static const char DATA_FILE_IDENT_[] = "dataFileIdent";
  constexpr static const char DATA_FILE_NAME_[] = "dataFileName";
  constexpr static const char DATA_PARSER_TYPE_[] = "dataFileType";
  // general purpose data file without type specific fields.
  constexpr static const char GENERAL_DATA_FILE_TYPE_[] = "generalFile";
  // VCF data file specific fields.
  constexpr static const char VCF_DATA_FILE_TYPE_[] = "vcfFile";
  constexpr static const char VCF_FILE_GENOME_[] =  "vcfGenome";
  constexpr static const char VCF_INFO_EVIDENCE_[] =  "vcfInfo";
  // VCF Info Evidence categories.
  constexpr static const char EVIDENCE_LIST_[] = "evidenceList";
  constexpr static const char EVIDENCE_IDENT_[] = "evidenceIdent";
  constexpr static const char EVIDENCE_INFO_LIST_[] = "vcfInfoList";
  constexpr static const char EVIDENCE_INFO_ITEM_[] = "vcfInfoItem";
   // Contig/Chromosome Alias categories.
  constexpr static const char ALIAS_LIST_[] = "aliasList";
  constexpr static const char ALIAS_IDENT_[] = "ident";
  constexpr static const char ALIAS_TYPE_[] = "chromosomeType";
  constexpr static const char ALIAS_ENTRY_[] = "alias";

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



#endif //KGL_PROPERTIES_H
