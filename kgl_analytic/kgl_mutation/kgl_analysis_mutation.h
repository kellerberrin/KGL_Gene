//
// Created by kellerberrin on 5/1/21.
//

#ifndef KGL_ANALYSIS_GENE_H
#define KGL_ANALYSIS_GENE_H


#include "kgl_analysis_virtual.h"
#include "kgl_Pf3k_COI.h"
#include "kgl_genealogy_parser.h"
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

  std::shared_ptr<const GenomeReference> ref_genome_ptr_;
  std::shared_ptr<const PopulationDB> population_ptr_;
  std::shared_ptr<const PopulationDB> unphased_population_ptr_;
  std::shared_ptr<const PopulationDB> clinvar_population_ptr_;
  std::shared_ptr<const GenomeAuxInfo> genome_aux_ptr_;
  std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr_;
  std::shared_ptr<const EnsemblHGNCResource> nomenclature_ptr_;
  // Results of the analysis. Type of gene membership is defined here.
  GenomeMutation gene_mutation_{VariantGeneMembership::BY_EXON};

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

  // From the OMIM entry #611162 available at https://www.omim.org/entry/611162
  inline static const std::vector<std::string> target_gene_map_ {
      "CD36", "GYPB", "FCGR2A", "FCGR2B", "ICAM1", "NCR3", "HBB", "NOS2", "TNF", "NCR3",
      "PROCR", "SLC4A1", "CISH", "MARVELD3", "FUT9", "VCAM1", "TIRAP", "GYPC", "FREM3", "GYPA",
      "G6PD", "CADM3",  "ACKR1", "ATP2B4", "CR1", "ABO", "HBA1", "BSG", "CD55", "LILRB1",
      "LAIR1" };

  // The list of genes to analyzed variant by variant. Must be ensembl codes (for now).
  inline static const std::vector<std::string> ensembl_gene_list_{

      "ENSG00000284690", "ENSG00000282992", "ENSG00000275019", "ENSG00000262576", "ENSG00000256797",
      "ENSG00000254521", "ENSG00000213402", "ENSG00000204345", "ENSG00000198178", "ENSG00000196371",
      "ENSG00000189184", "ENSG00000188211", "ENSG00000186407", "ENSG00000185475", "ENSG00000185187",
      "ENSG00000183840", "ENSG00000183019", "ENSG00000180549", "ENSG00000179213", "ENSG00000172794",
      "ENSG00000172322", "ENSG00000171840", "ENSG00000170956", "ENSG00000170425", "ENSG00000169704",
      "ENSG00000167850", "ENSG00000167123", "ENSG00000166589", "ENSG00000165682", "ENSG00000164713",
      "ENSG00000163600", "ENSG00000163485", "ENSG00000162897", "ENSG00000161649", "ENSG00000159674"

  };

};


} // namespace


#endif //KGL_ANALYSIS_GENE_H
