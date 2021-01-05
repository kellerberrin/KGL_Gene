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

  if (not getParameters(named_parameters, work_directory)) {

    ExecEnv::log().info("PfEMPAnalysis::initializeAnalysis; Analysis Id: {} problem parsing parameters", ident());
    return false;

  }

  reference_genomes_ = reference_genomes;
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



  std::shared_ptr<const LevenshteinLocal> levenshtein_distance_ptr(std::make_shared<const LevenshteinLocal>());
  std::shared_ptr<const Blosum80Local> blosum80_distance_ptr(std::make_shared<const Blosum80Local>());

  UPGMAMatrix upgma_matrix;

  VarGeneFamilyTree<AminoGeneDistance>( upgma_matrix,
                                        newick_file_name_,
                                        intron_file_name_,
                                        levenshtein_distance_ptr,
                                        reference_genomes_,
                                        "PFEMP1");

}



bool kgl::PfEMPAnalysis::getParameters(const ActiveParameterList& named_parameters, const std::string& work_directory) {


  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    if (block_vector.size() != 1) {

      ExecEnv::log().error("PfEMPAnalysis::getParameters; parameter block: {} vector size: {}, expected size = 1",
                           block_name, block_vector.size());
      return false;

    }

    ExecEnv::log().info("Analysis: {} parsing parameter block: {}", ident(), block_name);

    for (auto const& xml_vector : block_vector) {


      auto newick_opt = xml_vector.getString(NEWICK_FILE_);
      if (newick_opt) {

        newick_file_name_ = newick_opt.value().front();
        newick_file_name_ = Utility::filePath(newick_file_name_, work_directory);

      } else {

        ExecEnv::log().error("PfEMPAnalysis::getParameters; bad value for parameter: {}", NEWICK_FILE_);
        return false;

      }

      auto intron_opt = xml_vector.getString(INTRON_FILE_);
      if (intron_opt) {

        intron_file_name_ = intron_opt.value().front();
        intron_file_name_ = Utility::filePath(intron_file_name_, work_directory);

      } else {

        ExecEnv::log().error("PfEMPAnalysis::getParameters; bad value for parameter: {}", INTRON_FILE_);
        return false;

      }

    }

  }

  return true;

}

