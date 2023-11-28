//
// Created by kellerberrin on 3/1/21.
//

#ifndef KGL_ANALYSIS_PFEMP_H
#define KGL_ANALYSIS_PFEMP_H

#include "kgl_package_analysis_virtual.h"
#include "kgl_genome_collection.h"
#include "kgl_pf7_sample_parser.h"
#include "kgl_pf7_fws_parser.h"
#include "kga_analysis_lib_PfFilter.h"
#include "kgl_pf7_genetic_distance_parser.h"
#include "kgl_Pf7_physical_distance.h"
#include "kga_analysis_lib_seqmutation.h"
#include "kga_analysis_lib_seq_gene.h"

#include "kga_analysis_PfEMP_variant.h"
#include "kga_analysis_PfEMP_heterozygous.h"
#include "kga_analysis_PfEMP_FWS.h"

namespace kellerberrin::genome::analysis {   //  organization::project level namespace


class PfEMPAnalysis : public VirtualAnalysis {


public:

  PfEMPAnalysis() = default;

  ~PfEMPAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  constexpr static std::string IDENT {"PfEMP"};
  // Need a polymorphic version to interrogate VirtualAnalysis pointers.
  [[nodiscard]] std::string ident() const override { return IDENT; }
  // Simple creation factory function.
  [[nodiscard]] static std::unique_ptr<VirtualAnalysis> factory() { return std::make_unique<PfEMPAnalysis>(); }

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

  std::string ident_work_directory_;
  // Resources
  std::shared_ptr<const Pf7SampleResource> Pf7_sample_ptr_;
  std::shared_ptr<const Pf7FwsResource> Pf7_fws_ptr_;
  std::shared_ptr<const FilterPf7> filterPf7_ptr;
  std::shared_ptr<const Pf7GeneticDistanceResource> Pf7_genetic_distance_ptr_;
  constexpr static const double SAMPLE_LOCATION_RADIUS_{0.0};
  std::shared_ptr<const Pf7SampleLocation> Pf7_physical_distance_ptr_;
  constexpr static const std::string PF3D7_IDENT_{"Pf3D7_64"};
  std::shared_ptr<const GenomeReference> genome_3D7_ptr_;
  std::shared_ptr<const GenomeCollection> all_reference_genomes_ptr_;
  std::shared_ptr<MutateGenesReport> mutate_genes_ptr_; // Perform transcript level mutations for all genomes.


  // Gene family ident.
  constexpr static const std::string PFEMP1_FAMILY_{"PFEMP1"};
  constexpr static const std::string RUF6_FAMILY_{"RUF6"};
  constexpr static const std::string RIFIN_FAMILY_{"RIFIN"};
  constexpr static const std::string STEVOR_FAMILY_{"STEVOR"};
  constexpr static const std::string SURFIN_FAMILY_{"SURFIN"};
  constexpr static const std::string TRNA_FAMILY_{"TRNA"};
  constexpr static const std::string RIBOSOME_FAMILY_{"RIBOSOM"};


  // The genes we are interested in.
  GenomeGeneVariantAnalysis translation_gene_map_;   // tRNA, Ribosomes
  GenomeGeneVariantAnalysis antigenic_gene_map_;   // Var rifin stevor RUF6
  GenomeGeneVariantAnalysis all_gene_map_;     // All genes.

  // General variant statistics.
  HeteroHomoZygous hetero_homo_zygous_;
  CalcFWS calc_fws_;

  // Transcript Analysis
  AnalysisTranscriptFamily transcript_analysis_;
  constexpr static const std::string TRANSCRIPT_SUBDIRECTORY_{"transcript"};

  // File name constants.
  constexpr static const std::string VARIANT_COUNT_{"gene_variant_"};
  constexpr static const std::string VARIANT_COUNT_EXT_{".csv"};

  void checkDistanceMatrix(const std::shared_ptr<const PopulationDB> &all_population_ptr,
                           const std::shared_ptr<const PopulationDB> &filtered_population_ptr) const;

  void testPhysicalDistances();

  [[nodiscard]] GeneVector getAntiGenicGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);

  [[nodiscard]] GeneVector getAllGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);

  [[nodiscard]] GeneVector getTranslationGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);


};



} // namespace








#endif //KGL_KGL_ANALYSIS_PFEMP_H
