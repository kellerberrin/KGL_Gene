//
// Created by kellerberrin on 5/1/21.
//

#include "kgl_analysis_gene.h"
#include "kgl_analysis_gene_sequence.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::GeneAnalysis::initializeAnalysis(const std::string& work_directory,
                                              const ActiveParameterList& named_parameters,
                                              std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::GeneAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("Default VCF File Read for Analysis Id: {} called with Variant Population", ident());

  auto file_characteristic = data_ptr->dataCharacteristic();


  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::GeneAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::GeneAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}



void kgl::GeneAnalysis::performRegion() {

  std::string region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_561666";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   561666,
                                                   11753,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_462000";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_01_v3",
                                                   462000,
                                                   2000,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_0";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   0,
                                                   1200490,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_944950";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   944950,
                                                   135,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_941875";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   941875,
                                                   4293,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_569342";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   569342,
                                                   7467,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }

  region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_584668";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                                   "Pf3D7_04_v3",
                                                   584668,
                                                   7280,
                                                   population_ptr_,
                                                   genome_collection_ptr_->getGenome(analysis_genome),
                                                   region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }


}

