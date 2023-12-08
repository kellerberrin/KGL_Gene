//
// Created by kellerberrin on 7/11/20.
//

#include "kga_analysis_sequence.h"
#include "kga_analysis_lib_Pfgene.h"
#include "kga_analysis_lib_utility.h"
#include "kga_analysis_lib_PfFilter.h"

#include "kgl_variant_filter_db_variant.h"

#include <chrono>
#include <thread>

namespace kgl = kellerberrin::genome;
namespace kga = kellerberrin::genome::analysis;


// Setup the analytics to process VCF data.
bool kga::SequenceAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const ActiveParameterList& named_parameters,
                                               const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("SequenceAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }

  // Get subscribed resource for Pf3K COI data.
  Pf3KCOI_ptr_ = resource_ptr->getSingleResource<const Pf3kCOIResource>(ResourceProperties::PF3K_COI_RESOURCE_ID_);

  auto reference_genomes_ptr = std::make_shared<GenomeCollection>();
  for (auto const& genome_resource_ptr : resource_ptr->getResources(ResourceProperties::GENOME_RESOURCE_ID_)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());
    reference_genomes_ptr->addGenome(genome_ptr);

  }
  all_reference_genomes_ptr_ = reference_genomes_ptr;  // Assign to a pointer to const.

  auto pf3d7_opt = all_reference_genomes_ptr_->getOptionalGenome(PF3D7_IDENT_);
  if (not pf3d7_opt) {

    ExecEnv::log().critical("SequenceAnalysis::initializeAnalysis; Reference Genome: {} required for analysis - not supplied", PF3D7_IDENT_);

  }
  genome_3D7_ptr_ = pf3d7_opt.value();

  // Initialize the mutate object.
  mutate_genes_ptr_ = std::make_shared<MutateGenesReport>(genome_3D7_ptr_, SeqVariantFilterType::FRAMESHIFT_ADJUSTED, ident_work_directory_);

  // Do the UPGMA stuff.
  AnalysisGenePf::performGeneAnalysis(genome_3D7_ptr_, ident_work_directory_);

  // Analysis on RUf6 and PFEMP1.
//  transcript_analysis_.createAnalysisVector(KGAUtility::getRUF6Genes(genome_3D7_ptr_));
  transcript_analysis_.createAnalysisVector(KGAUtility::getPFEMP1Genes(genome_3D7_ptr_));

  return true;

}

// Perform the genetic analysis per iteration.
bool kga::SequenceAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> base_data_ptr) {

  ExecEnv::log().info("VCF File Read for Analysis Id: {} called with Variant Population", ident());

  // Superclass the population_ptr
  std::shared_ptr<const PopulationDB> population_ptr = std::dynamic_pointer_cast<const PopulationDB>(base_data_ptr);

  if (not population_ptr) {

    ExecEnv::log().error("Analysis: {}, expected a Population in file: {}", ident(), base_data_ptr->fileId());
    return false;

  }

  // Filter to COI = 1 only.
  FilterPf3k filter_coi(Pf3KCOI_ptr_);
  auto filtered_population_ptr = filter_coi.filterCOI(population_ptr);

  // Mutate all the relevant genes in the relevant contigs.
  mutate_genes_ptr_->mutatePopulation(filtered_population_ptr);

  transcript_analysis_.performFamilyAnalysis(filtered_population_ptr);

  return true;

}

// Perform the genetic analysis per iteration.
bool kga::SequenceAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kga::SequenceAnalysis::finalizeAnalysis() {

  mutate_genes_ptr_->printMutateReports();

  transcript_analysis_.printAllReports(ident_work_directory_, TRANSCRIPT_SUBDIRECTORY_);

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
