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

  select_gene_map_ = getSelectGeneMap(genome_3D7_ptr_);
  all_gene_map_ = getAllGeneMap(genome_3D7_ptr_);

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

  all_population_ptr_ = population_ptr;

  ExecEnv::log().info("Unfiltered Population: {}, Genome count: {}, Variant Count: {}",
                      all_population_ptr_->populationId(),
                      all_population_ptr_->getMap().size(),
                      all_population_ptr_->variantCount());

  filtered_population_ptr_ = Pf7_sample_ptr_->filterPassQCGenomes(all_population_ptr_);

  ExecEnv::log().info("QC Pass Filtered Population: {}, Genome count: {}, Variant Count: {}, Sample Data Count: {}",
                      filtered_population_ptr_->populationId(),
                      filtered_population_ptr_->getMap().size(),
                      filtered_population_ptr_->variantCount(),
                      Pf7_sample_ptr_->getMap().size());

  monoclonal_population_ptr_ = Pf7_fws_ptr_->filterFWS(FwsFilterType::GREATER_EQUAL, MONOCLONAL_FWS_THRESHOLD, filtered_population_ptr_);

  ExecEnv::log().info("MonoClonal Filtered Population: {}, Genome count: {}, Variant Count: {}, FWS Data Count: {}",
                      monoclonal_population_ptr_->populationId(),
                      monoclonal_population_ptr_->getMap().size(),
                      monoclonal_population_ptr_->variantCount(),
                      Pf7_fws_ptr_->getMap().size());

//  checkDistanceMatrix();

  getGeneVariants(select_gene_map_, filtered_population_ptr_);
  getGeneVariants(all_gene_map_, filtered_population_ptr_);

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
  writeGeneResults(select_gene_map_, variant_file_name);

  std::string all_variant_file_name = std::string(VARIANT_COUNT_) + "All" + std::string(VARIANT_COUNT_EXT_);
  all_variant_file_name = Utility::filePath(variant_file_name, ident_work_directory_);
  writeGeneResults(all_gene_map_, all_variant_file_name);

  return true;

}

void kgl::PfEMPAnalysis::writeGeneResults(const VariantGeneMap& gene_map, const std::string& variant_file_name) {

  std::ofstream variant_file(variant_file_name);

  if (not variant_file.good()) {

    ExecEnv::log().error("PfEMPAnalysis::writeGeneResults; Unable to open gene variant results file: {}", variant_file_name);
    return;

  }

  for (auto const& [gene_id, value_pair] : gene_map) {

    auto const& [gene_ptr, contig_ptr] = value_pair;

    auto sequence_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);

    size_t gene_length{0};
    if (not sequence_array_ptr->empty()) {

      gene_length = sequence_array_ptr->getFirst()->codingNucleotides();

    }

    double variant_rate{0.0};
    size_t variant_count = contig_ptr->variantCount();
    if (gene_length > 0) {

      variant_rate = static_cast<double>(variant_count) / static_cast<double>(gene_length);

    }


    variant_file << gene_ptr->contig()->contigId()
                 << CSV_DELIMITER_
                 << gene_id
                 << CSV_DELIMITER_
                 << gene_length
                 << CSV_DELIMITER_
                 << variant_count
                 << CSV_DELIMITER_
                 << variant_rate
                 << CSV_DELIMITER_
                 << gene_ptr->featureText(CSV_DELIMITER_)
                 << '\n';
  }

}

kgl::PfEMPAnalysis::VariantGeneMap kgl::PfEMPAnalysis::getSelectGeneMap(const std::shared_ptr<const GenomeReference>& genome_ptr) {

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

  // Place each of these genes into the GeneMap structure.
  VariantGeneMap gene_map;
  for (auto const& gene_ptr : var_gene_vector) {

    std::shared_ptr<ContigDB> gene_contig_ptr(std::make_shared<ContigDB>(gene_ptr->id()));

    std::pair<std::shared_ptr<const GeneFeature>, std::shared_ptr<ContigDB>> map_value{ gene_ptr, gene_contig_ptr};

    auto [iter, result] = gene_map.insert({gene_ptr->id(), map_value});

    if(not result) {

      ExecEnv::log().warn("PfEMPAnalysis::getSelectGeneMap; Duplicate Gene Id: {}", gene_ptr->id());

    }

  }

  return gene_map;

}



kgl::PfEMPAnalysis::VariantGeneMap kgl::PfEMPAnalysis::getAllGeneMap(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  VariantGeneMap gene_map;
  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [gene_id, gene_ptr] : contig_ptr->getGeneMap()) {

      std::shared_ptr<ContigDB> gene_contig_ptr(std::make_shared<ContigDB>(gene_ptr->id()));

      std::pair<std::shared_ptr<const GeneFeature>, std::shared_ptr<ContigDB>> map_value{ gene_ptr, gene_contig_ptr};

      auto [iter, result] = gene_map.insert({gene_ptr->id(), map_value});

      if(not result) {

        ExecEnv::log().warn("PfEMPAnalysis::getSelectGeneMap; Duplicate Gene Id: {}", gene_ptr->id());

      }

    }

  }

  return gene_map;

}


// Get variants only occurring within codingFeatures for all mRNA sequences.
void kgl::PfEMPAnalysis::getGeneVariants(VariantGeneMap& gene_map, const std::shared_ptr<const PopulationDB>& population_ptr) {

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      if (contig_ptr->getMap().empty()) {

        continue;

      }

      for (auto const &[gene_id, value_pair]: gene_map) {

        auto const &[gene_ptr, gene_contig_ptr] = value_pair;

        if (gene_ptr->contig()->contigId() == contig_id) {

          auto coding_sequence_array = GeneFeature::getTranscriptionSequences(gene_ptr);

          // Only analyze the first sequence.
          if (not coding_sequence_array->empty()) {

            auto sequence_ptr = coding_sequence_array->getFirst();

            for (const auto &[feature_id, feature_ptr]: sequence_ptr->getFeatureMap()) {

              gene_contig_ptr->merge(contig_ptr->subset(feature_ptr->sequence().begin(), feature_ptr->sequence().end()));

            } // For all cds

          } // if not empty

        } // if same contig.

      } // genemap

    } // contig

  } // genome

}

