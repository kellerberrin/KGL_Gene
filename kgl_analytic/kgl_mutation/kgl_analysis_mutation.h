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
#include "kgl_pubmed_api.h"



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
  // Various (low memory usage) requested resources
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
  // By Span is all variants with in intron+exon span of the gene
  // By Ensembl looks up the variants based on the vep ensembl code.
  // By Exon uses the gene exon addresses to find gene variants - warning assumes the first transcript.

  constexpr static const double FREQ_AFR_{0.0};
  constexpr static const double FREQ_ALL_{0.0};

  // Gene Allele Specific Analysis.
  GenerateGeneAllele gene_alleles_;
  // All Alleles With a Pubmed PMID identifier.
  GenerateGeneAllele all_pmid_alleles_;

  // Parameters and output files.
  std::string output_file_name_;
  std::string gene_allele_file_;
  std::string all_allele_file_;
  std::string literature_allele_file_;

  constexpr static const char* OUTPUT_FILE_ = "OutputFile";
  constexpr static const char* GENE_ALLELE_OUTPUT_FILE_ = "GeneAlleleOut";
  constexpr static const char* ALL_ALLELE_OUTPUT_FILE_ = "AllAlleleOut";
  constexpr static const char* LIT_ALLELE_FILE_ = "LitAllAllele";
  constexpr static const char OUTPUT_DELIMITER_ = ',';
  constexpr static const char* OUTPUT_FILE_EXT_ = ".csv";

  bool getParameters(const ActiveParameterList& named_parameters, const std::string& work_directory);

};


} // namespace


#endif //KGL_ANALYSIS_GENE_H
