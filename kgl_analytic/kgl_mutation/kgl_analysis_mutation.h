//
// Created by kellerberrin on 5/1/21.
//

#ifndef KGL_ANALYSIS_GENE_H
#define KGL_ANALYSIS_GENE_H


#include "kgl_analysis_virtual.h"
#include "kgl_Hsgenealogy_parser.h"
#include "kgl_uniprot_parser.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_analysis_mutation_gene_allele.h"


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

  // Parsed VCF Files.
  std::shared_ptr<const PopulationDB> population_ptr_;
  std::shared_ptr<const PopulationDB> unphased_population_ptr_;
  std::shared_ptr<const PopulationDB> clinvar_population_ptr_;
  // Various requested resources.
  std::shared_ptr<const GenomeReference> ref_genome_ptr_;
  std::shared_ptr<const HsGenomeAux> genome_aux_ptr_;
  std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr_;
  std::shared_ptr<const UniprotResource> uniprot_nomenclature_ptr_;
  std::shared_ptr<const EnsemblHGNCResource> ensembl_nomenclature_ptr_;
  // Results of the analysis. Type of gene membership is defined here.
  GenomeMutation gene_mutation_{VariantGeneMembership::BY_ENSEMBL};
  // By Span is all variants with in intron+exon span of the gene
  // By Ensembl looks up the variants based on the vep ensembl code.
  // By Exon uses the gene exon addresses to find gene variants - warning assumes the first transcript.

  constexpr static const double FREQ_AFR_{0.10};
  constexpr static const double UPPER_TAIL_AFR_{1E-05};   // 0.001 %
  constexpr static const double LOWER_TAIL_AFR_{1E-05};   // 0.001 %
  std::shared_ptr<GenerateGeneAllele> gene_allele_ptr_;

  std::string output_file_name_;
  std::string allele_output_file_name_;

  constexpr static const char* OUTPUT_FILE_ = "OutputFile";
  constexpr static const char* ALLELE_OUTPUT_FILE_ = "AlleleOut";
  constexpr static const char OUTPUT_DELIMITER_ = ',';
  constexpr static const char* OUTPUT_FILE_EXT_ = ".csv";

  bool getParameters(const ActiveParameterList& named_parameters, const std::string& work_directory);


};


} // namespace


#endif //KGL_ANALYSIS_GENE_H
