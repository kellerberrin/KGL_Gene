//
// Created by kellerberrin on 7/11/20.
//

#include "kgl_analysis_verify.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_variant_filter_db_offset.h"

#include <chrono>
#include <thread>

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::VerifyAnalysis::initializeAnalysis( const std::string& work_directory,
                                              const ActiveParameterList& named_parameters,
                                              const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident);

  }

  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("VerifyAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
                            ident_work_directory_);

  }

  auto reference_genomes_ptr = std::make_shared<GenomeCollection>();
  for (auto const& genome_resource_ptr : resource_ptr->getResources(ResourceProperties::GENOME_RESOURCE_ID_)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());
    reference_genomes_ptr->addGenome(genome_ptr);

  }
  all_reference_genomes_ptr_ = reference_genomes_ptr;  // Assign to a pointer to const.

  auto pf3d7_opt = all_reference_genomes_ptr_->getOptionalGenome(PF3D7_IDENT_);
  if (not pf3d7_opt) {

    ExecEnv::log().critical("VerifyAnalysis::initializeAnalysis; Reference Genome: {} required for analysis - not supplied", PF3D7_IDENT_);

  }
  genome_3D7_ptr_ = pf3d7_opt.value();

  // Initialize the mutate object.
  mutate_genes_ptr_ = std::make_shared<MutateGenes>(genome_3D7_ptr_);

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::VerifyAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> base_data_ptr) {

  ExecEnv::log().info("VCF File Read for Analysis Id: {} called with Variant Population", ident());

  // Superclass the population_ptr
  std::shared_ptr<const PopulationDB> population_ptr = std::dynamic_pointer_cast<const PopulationDB>(base_data_ptr);

  if (not population_ptr) {

    ExecEnv::log().error("Analysis: {}, expected a Population in file: {}", ident(), base_data_ptr->fileId());
    return false;

  }


  // Mutate all the relevant genes in the relevant contigs.
  ExecEnv::log().info("VerifyAnalysis::initializeAnalysis; Begin gene mutation");
  mutate_genes_ptr_->mutatePopulation(population_ptr);
  ExecEnv::log().info("VerifyAnalysis::initializeAnalysis; End gene mutation");

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::VerifyAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::VerifyAnalysis::finalizeAnalysis() {

  // Output the mutation statistics.
  std::string mutation_file_name = std::string("MutationTranscript") + std::string(VARIANT_COUNT_EXT_);
  mutation_file_name = Utility::filePath(mutation_file_name, ident_work_directory_);
  mutate_genes_ptr_->mutateAnalysis().printMutationTranscript(mutation_file_name);

  std::string validity_file_name = std::string("MutationValidity") + std::string(VARIANT_COUNT_EXT_);
  validity_file_name = Utility::filePath(validity_file_name, ident_work_directory_);
  mutate_genes_ptr_->mutateAnalysis().printMutationValidity(validity_file_name);

  mutation_file_name = std::string("MutationGenome") + std::string(VARIANT_COUNT_EXT_);
  mutation_file_name = Utility::filePath(mutation_file_name, ident_work_directory_);
  mutate_genes_ptr_->mutateAnalysis().printGenomeContig(mutation_file_name);

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
