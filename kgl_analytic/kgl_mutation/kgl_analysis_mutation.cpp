//
// Created by kellerberrin on 5/1/21.
//

#include "kgl_analysis_mutation.h"
#include "kgl_analysis_gene_sequence.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::MutationAnalysis::initializeAnalysis(const std::string& work_directory,
                                               const ActiveParameterList& named_parameters,
                                               const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter block: {}", ident(), parameter_ident);

  }

  auto genome_resource_vector = resource_ptr->getResources(RuntimeResourceType::GENOME_DATABASE);
  if (genome_resource_vector.size() != 1) {

    ExecEnv::log().critical("Analysis Id: {}, expected single (1) reference genome, actual count: {}", genome_resource_vector.size());

  }
  ref_genome_ptr_ = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_vector.front());

  auto ontology_resource_vector = resource_ptr->getResources(RuntimeResourceType::ONTOLOGY_DATABASE);
  if (ontology_resource_vector.size() != 1) {

    ExecEnv::log().critical("Analysis Id: {}, expected single (1) ontology database, actual count: {}", ontology_resource_vector.size());

  }

  ontology_db_ptr_ = std::dynamic_pointer_cast<const kol::OntologyDatabase>(ontology_resource_vector.front());

  auto nomenclature_resource_vector = resource_ptr->getResources(RuntimeResourceType::GENE_NOMENCLATURE);
  if (nomenclature_resource_vector.size() != 1) {

    ExecEnv::log().critical("Analysis Id: {}, expected single (1) Nomenclature database, actual count: {}", nomenclature_resource_vector.size());

  }

  nomenclature_ptr_ = std::dynamic_pointer_cast<const EnsemblHGNCResource>(nomenclature_resource_vector.front());

  if (not getParameters(named_parameters, work_directory)) {

    ExecEnv::log().critical("Analysis Id: {}, problem parsing parameters, program ends.", ident());

  }

  // Just in case.
  if (not ref_genome_ptr_ or not ontology_db_ptr_ or not nomenclature_ptr_) {

    ExecEnv::log().critical("Analysis Id: {}, A required resource is not defined, program ends.", ident());

  }

  return true;

}




bool kgl::MutationAnalysis::getParameters(const ActiveParameterList& named_parameters,
                                          const std::string& work_directory) {

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

  if (file_characteristic.data_structure == DataStructureEnum::PedGenome1000) {

    std::shared_ptr<const GenomePEDData> ped_data = std::dynamic_pointer_cast<const GenomePEDData>(data_ptr);

    if (ped_data) {

      ped_data_ = ped_data;
      ExecEnv::log().info("Analysis: {}, ped file: {} contains: {} PED records", ident(), ped_data->fileId(), ped_data->getMap().size());
      // Update the template populations.
      gene_mutation_.genomeAnalysis(ref_genome_ptr_, ped_data_, ontology_db_ptr_, nomenclature_ptr_);

      filterPedGenomes();

    } else {

      ExecEnv::log().critical("InbreedAnalysis::fileReadAnalysis, Analysis: {}, file: {} is not a PED Ancestor Object", ident(), data_ptr->fileId());
      return false;

    }

  }

  if (file_characteristic.data_implementation == DataImplEnum::PopulationVariant) {

    if (file_characteristic.data_source == DataSourceEnum::Genome1000
        or file_characteristic.data_source == DataSourceEnum::GnomadGenome3_1) {

      auto const_population = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

      if (not const_population) {

        ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast Genome1000 data file to population, severe error.");

      }

      auto population = std::const_pointer_cast<PopulationDB>(const_population);

      ExecEnv::log().info("Begin uniqueness filter for population: {} variant count: {}", population->populationId(), population->variantCount());
      auto pass_results = population->inSituFilter(PassFilter());
      auto diploid_results = population->inSituFilter(DiploidFilter());
      ExecEnv::log().info("Filtered Population: {} 'SNP and Pass' count: {}, 'Diploid' count: {}",
                          population->populationId(), pass_results.second, diploid_results.second);

      population_ptr_ = population;

      filterPedGenomes();


    } else if ( file_characteristic.data_source == DataSourceEnum::Gnomad2_1
               or file_characteristic.data_source == DataSourceEnum::Gnomad3_1) {

      unphased_population_ptr_ = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

      if (not unphased_population_ptr_) {

        ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast Gnomad 2.1 data file to population, severe error.");

      }

    } else if (file_characteristic.data_source == DataSourceEnum::Clinvar) {

      clinvar_population_ptr_ = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

      if (not clinvar_population_ptr_) {

        ExecEnv::log().critical("MutationAnalysis::fileReadAnalysis; Unable to cast Clinvar data file to population, severe error.");

      }

    }

  }

  return true;

}


void kgl::MutationAnalysis::filterPedGenomes() {


  if (population_ptr_ and ped_data_) {

    std::vector<GenomeId_t> ped_genomes;
    for (auto const& [genome, ped_data] :ped_data_->getMap()) {

      ped_genomes.push_back(genome);

    }

    size_t orginal_genomes = population_ptr_->getMap().size();
    auto non_const_population = std::const_pointer_cast<PopulationDB>(population_ptr_);
    auto count_pair = non_const_population->inSituFilter(GenomeFilter(ped_genomes));
    size_t final_genomes = non_const_population->getMap().size();
    ExecEnv::log().info("Filtered population: {} to Ped defined genomes, original count: {}, genomes: {}; filtered count: {}, genomes: {}",
                        population_ptr_->populationId(), count_pair.first, orginal_genomes, count_pair.second, final_genomes);

  }

}


// Perform the genetic analysis per iteration.
bool kgl::MutationAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());


  if (population_ptr_ and unphased_population_ptr_ and clinvar_population_ptr_ and ped_data_ and ontology_db_ptr_) {

    gene_mutation_.variantAnalysis(population_ptr_, unphased_population_ptr_, clinvar_population_ptr_, ped_data_);

  } else {

    ExecEnv::log().error("MutationAnalysis::iterationAnalysis; Necessary data files not defined for analysis");

  }

  std::pair<size_t, size_t> mem_pair = Utility::process_mem_usage2(); // pair.first is process vm_usage, pair.second is resident memory set.
  ExecEnv::log().info("Before Clear(), Variant Objects:{}, Data Blocks:{}, VM Usage: {}, Resident Memory: {}",
                      Variant::objectCount(), DataMemoryBlock::objectCount(), mem_pair.first, mem_pair.second);
  // Explicitly clean up the populations to recover memory.
  population_ptr_ = nullptr;
  unphased_population_ptr_ = nullptr;
  AuditMemory::trimFreeStore();
  mem_pair = Utility::process_mem_usage2();
  ExecEnv::log().info("After Clear(), Variant Objects:{}, Data Blocks:{}, VM Usage: {}, Resident Memory: {}",
                      Variant::objectCount(), DataMemoryBlock::objectCount(), mem_pair.first, mem_pair.second);

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::MutationAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  gene_mutation_.writeOutput(ped_data_, output_file_name_, OUTPUT_DELIMITER_);

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