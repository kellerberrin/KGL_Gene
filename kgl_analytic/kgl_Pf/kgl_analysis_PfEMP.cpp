//
// Created by kellerberrin on 3/1/21.
//

#include "kgl_analysis_PfEMP.h"
#include "kgl_upgma.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::PfEMPAnalysis::initializeAnalysis(const std::string& work_directory,
                                           const ActiveParameterList& named_parameters,
                                           std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  reference_genomes_ = reference_genomes;
  work_directory_ = work_directory;

  performPFEMP1UPGMA();

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB>) {

  ExecEnv::log().info("Default VCF File Read for Analysis Id: {} called with Variant Population", ident());

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::PfEMPAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}



void kgl::PfEMPAnalysis::performPFEMP1UPGMA() {


  std::string newick_file = Utility::filePath("newick_VAR", work_directory_) + ".txt";
  std::string intron_file = Utility::filePath("intron_VAR", work_directory_) + ".csv";

  std::shared_ptr<const LevenshteinLocal> levenshtein_distance_ptr(std::make_shared<const LevenshteinLocal>());
  std::shared_ptr<const Blosum80Local> blosum80_distance_ptr(std::make_shared<const Blosum80Local>());

  UPGMAMatrix upgma_matrix;

  VarGeneFamilyTree<AminoGeneDistance>( upgma_matrix,
                                        newick_file,
                                        intron_file,
                                        levenshtein_distance_ptr,
                                        reference_genomes_,
                                        "PFEMP1");

}



