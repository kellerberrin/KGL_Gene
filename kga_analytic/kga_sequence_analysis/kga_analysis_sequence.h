//
// Created by kellerberrin on 7/11/20.
//

#ifndef KGL_ANALYSIS_CHECK_H
#define KGL_ANALYSIS_CHECK_H


#include "kgl_package_analysis_virtual.h"
#include "kgl_genome_collection.h"
#include "kgl_pf3k_coi.h"
#include "kga_analysis_lib_seqmutation.h"
#include "kga_analysis_lib_seq_stats.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Analyzes Pf3k populations.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome::analysis {   //  organization::project level namespace


class SequenceAnalysis : public VirtualAnalysis {

public:

  SequenceAnalysis() = default;
  ~SequenceAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  constexpr static std::string IDENT {"PfSEQUENCE"};
  // Need a polymorphic version to interrogate VirtualAnalysis pointers.
  [[nodiscard]] std::string ident() const override { return IDENT; }
  // Simple creation factory function.
  [[nodiscard]] static std::unique_ptr<VirtualAnalysis> factory() { return std::make_unique<SequenceAnalysis>(); }

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

  constexpr static const std::string PF3D7_IDENT_{"Pf3D7_64"};
  std::shared_ptr<const GenomeReference> genome_3D7_ptr_;
  std::shared_ptr<const GenomeCollection> all_reference_genomes_ptr_;
  std::shared_ptr<MutateGenesReport> mutate_genes_ptr_; // Perform transcript level mutations for all genomes.

  // COI information for Pf3K samples.
  std::shared_ptr<const Pf3kCOIResource> Pf3KCOI_ptr_;

  // Transcript Analysis
  AnalysisTranscriptFamily transcript_analysis_;
  constexpr static const std::string TRANSCRIPT_SUBDIRECTORY_{"transcript"};
  constexpr static const std::string TRANSCRIPT_TREE_SUBDIRECTORY_{"transcript_tree"};

  constexpr static const std::string LABORATORY_3D7_TEXT_{"3D7_Laboratory"};  // 3D7 Lab. reference strain.
  constexpr static const std::string LABORATORY_DD2_TEXT_{"PG0008-CW"};  // DD2 Lab. reference strain.
  constexpr static const std::string LABORATORY_HB3_TEXT_{"PG0052-C"};  // HB3 Lab. reference strain.
  constexpr static const std::string LABORATORY_7G8_TEXT_{"7G8_7G8_Laborat"};  // 7GB Lab. reference strain.
  constexpr static const std::string LABORATORY_GB4_TEXT_{"GB4_GB4_Laborat"};  // GB4 Lab. reference strain.


  // constexpr static const std::string TREE_SELECTION_TEXT_{"Laboratory"};  // Only laboratory cultured strains.
  // constexpr static const std::string TREE_SELECTION_TEXT_{"Thailand"};
  constexpr static const std::string TREE_SELECTION_TEXT_{"Cambodia"};
  // constexpr static const std::string TREE_SELECTION_TEXT_{"Malawi"};

  // Allow for more complex selection criteria.
  constexpr static const auto sample_selection_ = [](const std::string& sample_text)->bool {

    return sample_text.contains(TREE_SELECTION_TEXT_)
    or sample_text.contains(LABORATORY_3D7_TEXT_)
    or sample_text.contains(LABORATORY_DD2_TEXT_)
    or sample_text.contains(LABORATORY_HB3_TEXT_)
    or sample_text.contains(LABORATORY_7G8_TEXT_)
    or sample_text.contains(LABORATORY_GB4_TEXT_);

  };

  std::string ident_work_directory_;
  constexpr static const std::string VARIANT_COUNT_EXT_{".csv"};

};



} // namespace


#endif //KGL_KGL_ANALYSIS_CHECK_H
