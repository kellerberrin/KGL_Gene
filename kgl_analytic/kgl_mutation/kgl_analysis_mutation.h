//
// Created by kellerberrin on 5/1/21.
//

#ifndef KGL_ANALYSIS_GENE_H
#define KGL_ANALYSIS_GENE_H


#include "kgl_analysis_virtual.h"
#include "kgl_hsgenealogy_parser.h"
#include "kgl_uniprot_parser.h"
#include "kgl_variant_sort_analysis.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_analysis_mutation_gene_allele.h"
#include "kgl_citation_parser.h"
#include "kgl_entrez_parser.h"
#include "kgl_bio_pmid_parser.h"
#include "kgl_pubmed_resource.h"



namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MutationAnalysis : public VirtualAnalysis {

public:

  MutationAnalysis() = default;
  ~MutationAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "MUTATION"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<MutationAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const ActiveParameterList& named_parameters,
                                         const std::shared_ptr<const AnalysisResources>& resource_ptr) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;

private:

  // These objects are cleaned up after each file (chromosome) iteration
  // to minimize memory use (which is still substantial).
  // Parsed VCF Files (reclaimed after each iteration).
  std::shared_ptr<const PopulationDB> population_ptr_;
  std::shared_ptr<const PopulationDB> unphased_population_ptr_;

  // All the objects below have the same lifetime as the analysis object.
  // These objects are considered (relatively) low memory usage.
  // Clinvar is relatively small, it is retained otherwise it would need to be
  // reloaded for each iteration.
  std::shared_ptr<const PopulationDB> clinvar_population_ptr_;
  // Various (relatively low memory usage) requested resources
  std::shared_ptr<const GenomeReference> ref_genome_ptr_;
  std::shared_ptr<const HsGenomeAux> genome_aux_ptr_;
  std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr_;
  std::shared_ptr<const UniprotResource> uniprot_nomenclature_ptr_;
  std::shared_ptr<const EnsemblHGNCResource> ensembl_nomenclature_ptr_;
  std::shared_ptr<const EntrezResource> entrez_nomenclature_ptr_;
  std::shared_ptr<const CitationResource> allele_citation_ptr_;
  std::shared_ptr<const PubmedRequester> pubmed_requestor_ptr_;

  // Results of the analysis. Type of gene membership is defined here.
  GenomeMutation gene_mutation_{VariantGeneMembership::BY_ENSEMBL};
  // By Span is all variants within the intron+exon span of the gene
  // By Ensembl looks up the variants based on the vep gene ensembl code.
  // By Exon uses the gene exon addresses to find gene variants - warning assumes the first gene transcript.

  // Gene Allele Specific Analysis.
  GenerateGeneAllele gene_alleles_;
  // All Alleles With a Pubmed PMID identifier.
  GenerateGeneAllele all_pmid_alleles_;
  // Profiles population variants with matching disease specific publications.
  GeneratePopulationAllele pop_pub_alleles_;

  // Parameters and output files.
  std::string ident_work_directory_;
  std::string population_lit_allele_directory_;
  std::string allele_literature_directory_;
  std::string allele_vep_directory_;

  const std::string ALLELE_LIT_SUBDIRECTORY_{"AlleleLiterature"};
  const std::string POPULATION_LIT_SUBDIRECTORY_{"PopulationLiterature"};
  const std::string ALLELE_VEP_SUBDIRECTORY_{"AlleleVEP"};

  std::string output_file_name_;

  constexpr static const char* OUTPUT_FILE_ = "OutputFile";
  constexpr static const char* GENE_ALLELE_OUTPUT_FILE_ = "_GeneAlleleOut";
  constexpr static const char* ALL_ALLELE_OUTPUT_FILE_ = "_AllAlleleOut";
  constexpr static const char* LIT_ALLELE_FILE_ = "_AllLitAllele";
  constexpr static const char* ALLELE_LIT_FILE_ = "_AllAlleleLit";
  constexpr static const char* POP_LIT_ALLELE_FILE_ = "_PopLitAllele";
  constexpr static const char OUTPUT_DELIMITER_ = ',';
  constexpr static const char* CSV_FILE_EXT_ = ".csv";
  constexpr static const char* TEXT_FILE_EXT_ = ".txt";

  bool getParameters(const ActiveParameterList& named_parameters);

};


} // namespace


#endif //KGL_ANALYSIS_GENE_H
