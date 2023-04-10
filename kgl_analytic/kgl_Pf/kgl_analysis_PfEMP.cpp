//
// Created by kellerberrin on 3/1/21.
//

#include "kgl_upgma.h"
#include "kgl_analysis_PfEMP.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::PfEMPAnalysis::initializeAnalysis(const std::string& work_directory,
                                            const ActiveParameterList& named_parameters,
                                            const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  Pf7_sample_ptr_ = resource_ptr->getSingleResource<const Pf7SampleResource>(ResourceProperties::PF7SAMPLE_RESOURCE_ID_);
  Pf7_fws_ptr_ = resource_ptr->getSingleResource<const Pf7FwsResource>(ResourceProperties::PF7FWS_RESOURCE_ID_);
  Pf7_distance_ptr_ = resource_ptr->getSingleResource<const Pf7DistanceResource>(ResourceProperties::PF7DISTANCE_RESOURCE_ID_);

  // Setup and clear the directories to hold analysis output.
  // The top level directory for this analysis type.
  ident_work_directory_ = work_directory + std::string("/") + ident();
  if (not Utility::createDirectory(ident_work_directory_)) {

    ExecEnv::log().critical("PfEMPAnalysis::initializeAnalysis, unable to create analysis results directory: {}",
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

    ExecEnv::log().critical("PfEMPAnalysis::initializeAnalysis; Reference Genome: {} required for analysis - not supplied", PF3D7_IDENT_);

  }
  genome_3D7_ptr_ = pf3d7_opt.value();

//  performPFEMP1UPGMA();

  // Select the genes we are interested in analyzing for genome variants.
  select_gene_map_.setGeneVector(getSelectGenes(genome_3D7_ptr_));
  all_gene_map_.setGeneVector(getAllGenes(genome_3D7_ptr_));

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::PfEMPAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> base_data_ptr) {

  ExecEnv::log().info("File Read for Analysis Id: {} called for file: {}", ident(), base_data_ptr->fileId());

  // Superclass the population_ptr
  std::shared_ptr<const PopulationDB> population_ptr = std::dynamic_pointer_cast<const PopulationDB>(base_data_ptr);

  if (not population_ptr) {

    ExecEnv::log().error("Analysis: {}, expected a Population in file: {}", ident(), base_data_ptr->fileId());
    return false;

  }

  ExecEnv::log().info("Unfiltered Population: {}, Genome count: {}, Variant Count: {}",
                      population_ptr->populationId(),
                      population_ptr->getMap().size(),
                      population_ptr->variantCount());

  auto filtered_population_ptr = Pf7_sample_ptr_->filterPassQCGenomes(population_ptr);

  ExecEnv::log().info("QC Pass Filtered Population: {}, Genome count: {}, Variant Count: {}, Sample Data Count: {}",
                      filtered_population_ptr->populationId(),
                      filtered_population_ptr->getMap().size(),
                      filtered_population_ptr->variantCount(),
                      Pf7_sample_ptr_->getMap().size());

  auto monoclonal_population_ptr = Pf7_fws_ptr_->filterFWS(FwsFilterType::GREATER_EQUAL, MONOCLONAL_FWS_THRESHOLD, filtered_population_ptr);

  ExecEnv::log().info("MonoClonal Filtered Population: {}, Genome count: {}, Variant Count: {}, FWS Data Count: {}",
                      monoclonal_population_ptr->populationId(),
                      monoclonal_population_ptr->getMap().size(),
                      monoclonal_population_ptr->variantCount(),
                      Pf7_fws_ptr_->getMap().size());

//  checkDistanceMatrix(population_ptr, filtered_population_ptr);

  select_gene_map_.getGeneVariants(filtered_population_ptr);
  all_gene_map_.getGeneVariants(filtered_population_ptr);

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

  std::string variant_file_name = std::string(VARIANT_COUNT_) + "Select" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  select_gene_map_.writeGeneResults(variant_file_name);

  std::string all_variant_file_name = std::string(VARIANT_COUNT_) + "All" + std::string(VARIANT_COUNT_EXT_);
  all_variant_file_name = Utility::filePath(all_variant_file_name, ident_work_directory_);
  all_gene_map_.writeGeneResults(all_variant_file_name);

  return true;

}



kgl::GeneVector kgl::PfEMPAnalysis::getSelectGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  // Get the gene families of interest.
  auto var_gene_vector = getGeneVector(genome_ptr, PFEMP1_FAMILY_);
  auto ruf6_gene_vector = getGeneVector(genome_ptr, RUF6_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), ruf6_gene_vector.begin(), ruf6_gene_vector.end() );
  auto rifin_gene_vector = getGeneVector(genome_ptr, RIFIN_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), rifin_gene_vector.begin(), rifin_gene_vector.end() );
  auto stevor_gene_vector = getGeneVector(genome_ptr, STEVOR_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), stevor_gene_vector.begin(), stevor_gene_vector.end() );
  auto surfin_gene_vector = getGeneVector(genome_ptr, SURFIN_FAMILY_);
  var_gene_vector.insert(var_gene_vector.end(), surfin_gene_vector.begin(), surfin_gene_vector.end() );
  //  auto ncRNA_gene_vector = getncRNAGeneVector(genome_ptr);
  //  var_gene_vector.insert( var_gene_vector.end(), ncRNA_gene_vector.begin(), ncRNA_gene_vector.end() );

  return var_gene_vector;

}


kgl::GeneVector kgl::PfEMPAnalysis::getAllGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  GeneVector all_genes;
  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [gene_id, gene_ptr] : contig_ptr->getGeneMap()) {

      all_genes.push_back(gene_ptr);

    }

  }

  return all_genes;

}

