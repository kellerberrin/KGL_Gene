//
// Created by kellerberrin on 3/1/21.
//

#ifndef KGL_ANALYSIS_PFEMP_H
#define KGL_ANALYSIS_PFEMP_H


#include "kgl_analysis_virtual.h"
#include "kgl_genome_collection.h"
#include "kgl_pf7_sample_parser.h"
#include "kgl_pf7_fws_parser.h"
#include "kgl_pf7_genetic_distance_parser.h"
#include "kgl_Pf7_physical_distance.h"
#include "kgl_analysis_PfEMP_variant.h"
#include "kgl_analysis_PfEMP_heterozygous.h"
#include "kgl_analysis_PfEMP_FWS.h"
#include "kgl_variant_filter_features.h"
#include "kgl_analysis_PfEMP_overlap.h"
#include "kgl_analysis_PfEMP_mutation.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class PfEMPAnalysis : public VirtualAnalysis {


public:

  PfEMPAnalysis() = default;

  ~PfEMPAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return PFEMPANALYSIS_IDENT_; }

  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<PfEMPAnalysis>(); }

  // Setup the analytics to process VCF data.
  [[nodiscard]] bool initializeAnalysis(const std::string &work_directory,
                                        const ActiveParameterList &named_parameters,
                                        const std::shared_ptr<const AnalysisResources> &resource_ptr) override;


  // Perform the genetic analysis per VCF file
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] bool finalizeAnalysis() override;


private:

  constexpr static const char PFEMPANALYSIS_IDENT_[]{"PfEMP"};

  std::string ident_work_directory_;
  // Resources
  std::shared_ptr<const Pf7SampleResource> Pf7_sample_ptr_;
  std::shared_ptr<const Pf7FwsResource> Pf7_fws_ptr_;
  std::shared_ptr<const Pf7GeneticDistanceResource> Pf7_genetic_distance_ptr_;
  constexpr static const double SAMPLE_LOCATION_RADIUS_{0.0};
  std::shared_ptr<const Pf7SampleLocation> Pf7_physical_distance_ptr_;
  constexpr static const char PF3D7_IDENT_[]{"Pf3D7_62"};
  std::shared_ptr<const GenomeReference> genome_3D7_ptr_;
  std::shared_ptr<const GenomeCollection> all_reference_genomes_ptr_;
  std::shared_ptr<const MutateGenes> mutate_genes_ptr_;

  // Filter constants.
  constexpr static const bool CODING_FILTER_ACTIVE_{true}; // Restrict to gene coding areas only.
  constexpr static const bool FILTER_QC_ACTIVE_{true};
  constexpr static const bool FILTER_FWS_ACTIVE_{false};
  constexpr static const bool MLEAF_FILTER_ACTIVE_{false};
  constexpr static const bool AF_FILTER_ACTIVE_{false};
  constexpr static const bool SNP_FILTER_ACTIVE_{true};  // Only SNP variants
  constexpr static const double VARIANT_FREQUENCY_CUTOFF_{0.0};

  void performPFEMP1UPGMA();

  // Available Parameters
  constexpr static const char *NEWICK_FILE_ = "NewickFile";
  constexpr static const char *INTRON_FILE_ = "IntronFile";

  // Intron promoter sequences.
  constexpr static const char *I_PROMOTER_ = "TGTATGTG";
  constexpr static const char *I_COMPLEMENT_PROMOTER_ = "ACATACAC";
  constexpr static const char *I_5_PROMOTER_ = "TCATA";

  // Min seq size for UPGMA analysis.
  constexpr static const size_t MIN_SEQUENCE_LENGTH_ = 10;
  // Gene family ident.
  constexpr static const char *PFEMP1_FAMILY_ = "PFEMP1";
  constexpr static const char *RUF6_FAMILY_ = "RUF6";
  constexpr static const char *RIFIN_FAMILY_ = "RIFIN";
  constexpr static const char *STEVOR_FAMILY_ = "STEVOR";
  constexpr static const char *SURFIN_FAMILY_ = "SURFIN";
  constexpr static const char *TRNA_FAMILY_ = "TRNA";
  constexpr static const char *RIBOSOME_FAMILY_ = "RIBOSOM";


  // The genes we are interested in.
  GenomeGeneVariantAnalysis translation_gene_map_;   // tRNA, Ribosomes
  GenomeGeneVariantAnalysis antigenic_gene_map_;   // Var rifin stevor RUF6
  GenomeGeneVariantAnalysis all_gene_map_;     // All genes.

  // General variant statistics.
  HeteroHomoZygous hetero_homo_zygous_;
  CalcFWS calc_fws_;

  // Check overlapping genes.
  std::shared_ptr<OverlapGenes> overlap_ptr_;

  // File name constants.
  constexpr static const char *NEWICK_{"newick_"};
  constexpr static const char *NEWICK_EXT_{".txt"};
  constexpr static const char *INTRON_{"intron_"};
  constexpr static const char *INTRON_EXT_{".csv"};
  constexpr static const char *VARIANT_COUNT_{"gene_variant_"};
  constexpr static const char *VARIANT_COUNT_EXT_{".csv"};
  constexpr static const char CSV_DELIMITER_ = ',';

  // Return a vector of genes that have a particular text fragment in the description.
  [[nodiscard]] GeneVector getGeneVector(const std::shared_ptr<const GenomeReference> &genome_ptr, const std::string &description_text) const;

  [[nodiscard]] GeneVector getncRNAGeneVector(const std::shared_ptr<const GenomeReference> &genome_ptr,
                                              const std::string &desc_uc_text = "",
                                              size_t max_size = 0) const;

  [[nodiscard]] GeneVector proximityGenes(size_t radius,
                                          const std::shared_ptr<const GeneFeature> &target_ptr,
                                          const GeneVector &gene_vector) const;

  // Analyze the introns of Var family genes.
  void varIntron(const GeneVector &gene_vector, const std::string &intron_file) const;

  void geneFamilyUPGMA(const std::shared_ptr<const GenomeReference> &genome_ptr,
                       const GeneVector &gene_vector,
                       const std::string &upgma_file_name,
                       const std::string &family_text) const;

  void checkDistanceMatrix(const std::shared_ptr<const PopulationDB> &all_population_ptr,
                           const std::shared_ptr<const PopulationDB> &filtered_population_ptr) const;

  void testPhysicalDistances();

  [[nodiscard]] GeneVector getAntiGenicGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);

  [[nodiscard]] GeneVector getAllGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);

  [[nodiscard]] GeneVector getTranslationGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);

  // Quality filter the variants using read depth, VQSLOD and other statistics
  std::shared_ptr<kgl::PopulationDB> qualityFilter(const std::shared_ptr<const PopulationDB> &unfiltered_population);

  // Perform mutation analysis of ncRNA and Protein genes.
  void performMutation(const std::shared_ptr<const PopulationDB> &filtered_population);

};



} // namespace








#endif //KGL_KGL_ANALYSIS_PFEMP_H
