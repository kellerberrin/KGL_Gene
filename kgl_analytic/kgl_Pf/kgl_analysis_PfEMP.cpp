//
// Created by kellerberrin on 3/1/21.
//

#include "kgl_analysis_PfEMP.h"
#include "kgl_variant_filter_features.h"

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
  Pf7_genetic_distance_ptr_ = resource_ptr->getSingleResource<const Pf7GeneticDistanceResource>(ResourceProperties::PF7DISTANCE_RESOURCE_ID_);
  Pf7_physical_distance_ptr_ = std::make_shared<const Pf7SampleLocation>(*Pf7_sample_ptr_);

  testPhysicalDistances();

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

  // Initialize the overlap analysis.
  overlap_ptr_ = std::make_shared<OverlapGenes>(genome_3D7_ptr_);

  // Initialize the mutate object.
  mutate_genes_ptr_ = std::make_shared<const MutateGenes>(genome_3D7_ptr_);

//  performPFEMP1UPGMA();

  // Select the genes we are interested in analyzing for genome variants.
  antigenic_gene_map_.setGeneVector(getAntiGenicGenes(genome_3D7_ptr_));
  all_gene_map_.setGeneVector(getAllGenes(genome_3D7_ptr_));
  translation_gene_map_.setGeneVector(getTranslationGenes(genome_3D7_ptr_));

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

  auto filtered_population_ptr = qualityFilter(population_ptr);

  ExecEnv::log().info("Population Returned Filtered Size Genome count: {}, Variant Count: {}",
                      filtered_population_ptr->getMap().size(),
                      filtered_population_ptr->variantCount());

//  checkDistanceMatrix(population_ptr, filtered_population_ptr);

  // Analyze gene variant info
  translation_gene_map_.getGeneVariants(filtered_population_ptr);
  antigenic_gene_map_.getGeneVariants(filtered_population_ptr);
  all_gene_map_.getGeneVariants(filtered_population_ptr);

  // Analyze for Homozygous and overlapping variants.
  hetero_homo_zygous_.analyzeVariantPopulation(filtered_population_ptr, Pf7_fws_ptr_, Pf7_sample_ptr_);

  // Calculate the FWS statistics.
  calc_fws_.calcFwsStatistics(filtered_population_ptr);

  // Mutate all the relevant genes in the relevant contigs.
  performMutation(filtered_population_ptr);

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

  // Print results files.

  std::string variant_file_name = std::string(VARIANT_COUNT_) + "Translation" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  translation_gene_map_.writeGeneResults(variant_file_name);

  variant_file_name = std::string(VARIANT_COUNT_) + "AntiGenic" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  antigenic_gene_map_.writeGeneResults(variant_file_name);

  variant_file_name = std::string(VARIANT_COUNT_) + "All" + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  all_gene_map_.writeGeneResults(variant_file_name);


  // Get location summaries
  auto location_summary_map = hetero_homo_zygous_.location_summary(Pf7_sample_ptr_,
                                                                   Pf7_physical_distance_ptr_,
                                                                   SAMPLE_LOCATION_RADIUS_,
                                                                   Pf7_fws_ptr_);

  // Update the FIS statistic.
  hetero_homo_zygous_.UpdateSampleLocation(location_summary_map);

  variant_file_name = std::string("VariantStatistics") + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  hetero_homo_zygous_.write_variant_results(variant_file_name, location_summary_map);

  variant_file_name = std::string("VariantLocation") + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);

  hetero_homo_zygous_.write_location_results( variant_file_name, location_summary_map);

  variant_file_name = std::string("GenomeFWS") + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  calc_fws_.writeGenomeResults(Pf7_fws_ptr_, variant_file_name);

  variant_file_name = std::string("VariantFWS") + std::string(VARIANT_COUNT_EXT_);
  variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  calc_fws_.writeVariantResults(variant_file_name);

  // Output the mutation statistics.
  std::string mutation_file_name = std::string("MutationTranscript") + std::string(VARIANT_COUNT_EXT_);
  mutation_file_name = Utility::filePath(mutation_file_name, ident_work_directory_);
  mutate_analysis_.printMutationTranscript(mutation_file_name);
  mutation_file_name = std::string("MutationGenome") + std::string(VARIANT_COUNT_EXT_);
  mutation_file_name = Utility::filePath(mutation_file_name, ident_work_directory_);
  mutate_analysis_.printGenomeContig(mutation_file_name);

  // Overlap to log file.
  overlap_ptr_->printResults();

  return true;

}



kgl::GeneVector kgl::PfEMPAnalysis::getAntiGenicGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

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

kgl::GeneVector kgl::PfEMPAnalysis::getTranslationGenes(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  // Get the gene families of interest.
  auto trna_gene_vector = getncRNAGeneVector(genome_ptr, TRNA_FAMILY_);
  auto ribosome_gene_vector = getGeneVector(genome_ptr, RIBOSOME_FAMILY_);
  trna_gene_vector.insert(trna_gene_vector.end(), ribosome_gene_vector.begin(), ribosome_gene_vector.end() );

  return trna_gene_vector;

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

void kgl::PfEMPAnalysis::testPhysicalDistances() {

  ExecEnv::log().info("Number of Locations: {}", Pf7_physical_distance_ptr_->locationMap().size());

  const double radius{1000.0};
  const std::string check_location{"Greater Accra"};
  for (auto const& [location, location_record] : Pf7_physical_distance_ptr_->locationMap()) {

    auto const& [loc, type] = location_record.location();
    ExecEnv::log().info("Location: {}, Location type: {}, Samples: {}, latitude: {}, longitude: {}, Distance from {}: {}, Radius: {}, Locations:{}, Samples: {}",
                        location, (type == LocationType::City ? "City" : "Country"),
                        location_record.locationSamples().size(),
                        location_record.latitudeDegrees(), location_record.longitudeDegrees(),
                        check_location, Pf7_physical_distance_ptr_->calculateDistance(check_location, location),
                        radius, Pf7_physical_distance_ptr_->locationRadius(location, radius).size(),
                        Pf7_physical_distance_ptr_->sampleRadius(location, radius).size());

  }

  const double radius2{0.0};
  const std::string check_country{"Ghana"};
  for (auto const& [location, location_record] : Pf7_physical_distance_ptr_->locationMap()) {

    auto const& [loc, type] = location_record.location();
    ExecEnv::log().info("Location: {}, Location type: {}, Samples: {}, latitude: {}, longitude: {}, Distance from {}: {}, Radius: {}, Locations:{}, Samples: {}",
                        location, (type == LocationType::City ? "City" : "Country"),
                        location_record.locationSamples().size(),
                        location_record.latitudeDegrees(), location_record.longitudeDegrees(),
                        check_country, Pf7_physical_distance_ptr_->distance(check_country, location),
                        radius2, Pf7_physical_distance_ptr_->locationRadius(location, radius2).size(),
                        Pf7_physical_distance_ptr_->sampleRadius(location, radius2).size());

  }



}

void kgl::PfEMPAnalysis::performMutation(const std::shared_ptr<const PopulationDB> &filtered_population_ptr) {

  ExecEnv::log().info("PfEMPAnalysis::performMutation; Begin gene mutation");

  // Get the active contigs in this population.
  auto contig_map = filtered_population_ptr->contigCount();

  for (auto const& [contig_id, variant_count] : contig_map) {

    if (variant_count > 0) {

      auto gene_vector = mutate_genes_ptr_->contigGenes(contig_id);
      ExecEnv::log().info("PfEMPAnalysis::performMutation; Mutating contig: {}, gene count: {}", contig_id, gene_vector.size());

      for (auto const& gene_ptr : gene_vector) {

        auto transcription_array = GeneFeature::getTranscriptionSequences(gene_ptr);
        for (auto const& [transcript_id,  transcript_ptr] : transcription_array->getMap()) {

          auto const [total_variants, duplicate_variants] = mutate_genes_ptr_->mutateTranscript( gene_ptr,
                                                                                                 transcript_id,
                                                                                                 filtered_population_ptr,
                                                                                                 genome_3D7_ptr_);

          mutate_analysis_.addTranscriptRecord(TranscriptMutateRecord(gene_ptr, transcript_ptr, total_variants, duplicate_variants));

        }

      }

    }

  }

  ExecEnv::log().info("PfEMPAnalysis::performMutation; End gene mutation");

}

