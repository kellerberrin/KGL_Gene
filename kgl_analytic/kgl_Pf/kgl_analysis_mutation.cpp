//
// Created by kellerberrin on 5/1/21.
//

#include "kgl_analysis_mutation.h"
#include "kgl_analysis_gene_sequence.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const ActiveParameterList& named_parameters,
                                               std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  if (reference_genomes->getMap().size() != 1) {

    ExecEnv::log().critical("Analysis Id: {}, expected single (1) reference genome, actual count: {}", reference_genomes->getMap().size());

  }

  auto [genome_id, genome_ptr] = *(reference_genomes->getMap().begin());
  ref_genome_ptr_ = genome_ptr;

  if (not getParameters(named_parameters, work_directory)) {

    ExecEnv::log().critical("Analysis Id: {}, problem parsing parameters, program ends.", ident());

  }

  // Sanity Check.
  if (not ref_genome_ptr_) {

    ExecEnv::log().critical("Analysis Id: {}, no genomes defined, program ends.", ident());

  }

  gene_mutation_.genomeAnalysis(ref_genome_ptr_);

  return true;

}




bool kgl::MutationAnalysis::getParameters(const ActiveParameterList& named_parameters,
                                          const std::string& work_directory) {

  work_directory_ = work_directory;

  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    if (block_vector.size() != 1) {

      ExecEnv::log().error("MutationAnalysis::getParameters; parameter block: {} vector size: {}, expected size = 1",
                           block_name, block_vector.size());
      return false;

    }

    ExecEnv::log().info("Analysis: {} parsing parameter block: {}", ident(), block_name);

    for (auto const& xml_vector : block_vector) {


      auto output_opt = xml_vector.getString(OUTPUT_FILE_);
      if (output_opt) {

        output_file_name_ = output_opt.value().front() + std::string(OUTPUT_FILE_EXT_);
        output_file_name_ = Utility::filePath(output_file_name_, work_directory);
        ExecEnv::log().info("Analysis: {} outputfile: {}", ident(), output_file_name_);

      } else {

        ExecEnv::log().error("MutationAnalysis::getParameters; bad value for parameter: {}", OUTPUT_FILE_);
        return false;

      }

    }

  }

  return true;

}



// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called with file: {}", ident(), data_ptr->fileId());

  auto file_characteristic = data_ptr->dataCharacteristic();

  if (file_characteristic.data_source == DataSourceEnum::Pf3kCOI) {

    pf3k_coi_ptr_ = std::dynamic_pointer_cast<const Pf3kCOIDB>(data_ptr);

    if (pf3k_coi_ptr_) {

      ExecEnv::log().info("Processed Pf3k complexity of infection file: {}", pf3k_coi_ptr_->fileId());

    } else {

      ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast to the Pf3K complexity of infection file, severe error.");

    }

  }

  if (file_characteristic.data_structure == DataStructureEnum::PedGenome1000) {

    std::shared_ptr<const GenomePEDData> ped_data = std::dynamic_pointer_cast<const GenomePEDData>(data_ptr);

    if (ped_data) {

      ped_data_ = ped_data;
      ExecEnv::log().info("Analysis: {}, ped file: {} contains: {} PED records", ident(), ped_data->fileId(), ped_data->getMap().size());

    } else {

      ExecEnv::log().critical("InbreedAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a PED Ancestor Object", ident(), data_ptr->fileId());
      return false;

    }

  }

  if (file_characteristic.data_implementation == DataImplEnum::PopulationVariant) {

    population_ptr_ = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

    if (not population_ptr_) {

      ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast data file to population, severe error.");

    }

  }

  return true;

}


std::shared_ptr<const kgl::PopulationDB> kgl::MutationAnalysis::createUnphased(const std::shared_ptr<const PopulationDB>& population) {

  ExecEnv::log().info("InbreedAnalysis::processDiploid; Creating unique unphased population using Diploid Population.");
  std::shared_ptr<GenomeDB> unphased_genome_ptr = population->uniqueUnphasedGenome();
  std::shared_ptr<PopulationDB> unphased_unique_ptr = std::make_shared<PopulationDB>(population->populationId(), population->dataSource());
  unphased_unique_ptr->addGenome(unphased_genome_ptr);
  ExecEnv::log().info("InbreedAnalysis::processDiploid; Created unique unphased population, variant count: {}.", unphased_unique_ptr->variantCount());

  return unphased_unique_ptr;

}


// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());


  if (population_ptr_) {

    auto file_characteristic = population_ptr_->dataCharacteristic();

    if (file_characteristic.data_source == DataSourceEnum::Genome1000 and ped_data_) {

      gene_mutation_.variantAnalysis100(population_ptr_, ped_data_);

    } else {

      gene_mutation_.variantAnalysis(population_ptr_);

    }

  } else {

    ExecEnv::log().error("MutationAnalysis::iterationAnalysis; No deta files defined for analysis");

  }

  population_ptr_ = nullptr;

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  gene_mutation_.writeOutput100(output_file_name_, OUTPUT_DELIMITER_);

  return true;

}

/*

void kgl::MutationAnalysis::performRegion() {

  const GenomeId_t analysis_genome = "Pf3D7_47";

  std::string region_fasta_file = "malawi_fb_SRR609075";
  region_fasta_file += "_561666";
  region_fasta_file = Utility::filePath(region_fasta_file, work_directory_) + ".fasta";
  if (not GenomicSequence::mutateGenomeRegion("malawi_fb_SRR609075",
                                              "Pf3D7_04_v3",
                                              561666,
                                              11753,
                                              population_ptr_,
                                              ref_genome_ptr_->getGenome(analysis_genome),
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
                                              ref_genome_ptr_->getGenome(analysis_genome),
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
                                              ref_genome_ptr_->getGenome(analysis_genome),
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
                                              ref_genome_ptr_->getGenome(analysis_genome),
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
                                              ref_genome_ptr_->getGenome(analysis_genome),
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
                                              ref_genome_ptr_->getGenome(analysis_genome),
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
                                              ref_genome_ptr_->getGenome(analysis_genome),
                                              region_fasta_file)) {

    ExecEnv::log().error("PhylogeneticAnalysis::performRegion(); analysis fails");

  }


}

*/