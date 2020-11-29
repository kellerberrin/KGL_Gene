//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_inbreed_locus.h"
#include "kgl_analysis_mutation_synthetic.h"
#include "kgl_filter.h"


#include <fstream>

namespace kgl = kellerberrin::genome;



bool kgl::InbreedingAnalysis::InbreedingAll( std::shared_ptr<const GenomeReference> genome_ptr,
                                             std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                             std::shared_ptr<const DiploidPopulation> diploid_ptr,
                                             std::shared_ptr<const GenomePEDData> ped_data_ptr,
                                             const std::string& output_file_name,
                                             InbreedingParameters& parameters) {

  static const std::string gnomad_3_1_fragment = "r3.1";
  // Setup the runtime parameters.

  parameters.outputFile(output_file_name);
  std::string file_name =  output_file_name + "_" + parameters.inbreedingAlgorthim() + "_";
  parameters.outputFile(file_name);

  // Only the diploid population presented.
  // Create a synthetic unphased population using the diploid population.
  if (diploid_ptr and ped_data_ptr and not unphased_ptr) {

    ExecEnv::log().info("InbreedingAnalysis::Inbreeding; Creating unique unphased population using 1000 Genomes.");
    std::shared_ptr<GenomeVariant> unphased_genome_ptr = diploid_ptr->uniqueUnphasedGenome<GenomeVariant>();
    std::shared_ptr<UnphasedPopulation> unphased_unique_ptr = std::make_shared<UnphasedPopulation>("UniqueUnphasedPopulation");
    unphased_unique_ptr->addGenome(unphased_genome_ptr);
    ExecEnv::log().info("InbreedingAnalysis::Inbreeding; Created unique unphased population, variant count: {}.", unphased_unique_ptr->variantCount());
    unphased_ptr = unphased_unique_ptr;
    parameters.lociiArguments().variantSource(VariantDatabaseSource::GENOMES_1000);

  } else if (unphased_ptr) {

    // Check if Gnomad 3.1 else Gnomad 2.1
    size_t find_pos = unphased_ptr->populationId().find(gnomad_3_1_fragment);
    if (find_pos != std::string::npos) {

      parameters.lociiArguments().variantSource(VariantDatabaseSource::GNOMAD3_1);

    } else {

      parameters.lociiArguments().variantSource(VariantDatabaseSource::GNOMAD2_1);

    }

  } else {

    ExecEnv::log().error("MutationAnalysis::iterationAnalysis, missing required data files");
    return false;

  }

  if (not InbreedingAnalysis::Inbreeding( genome_ptr,
                                          unphased_ptr,
                                          diploid_ptr,
                                          ped_data_ptr,
                                          parameters)) {

    ExecEnv::log().error("MutationAnalysis::iterationAnalysis, problem with population inbreeding analysis, algorithm: {}",
                         parameters.inbreedingAlgorthim());
    return false;

  }


  return true;

}

bool kgl::InbreedingAnalysis::Inbreeding( std::shared_ptr<const GenomeReference> genome_ptr,
                                          std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                          std::shared_ptr<const DiploidPopulation> diploid_ptr,
                                          std::shared_ptr<const GenomePEDData> ped_data_ptr,
                                          InbreedingParameters& parameters) {


  if (diploid_ptr and unphased_ptr and ped_data_ptr) {

    ExecEnv::log().info("Filtered Population : {}  and Joined Population: {}  both active",
                        diploid_ptr->populationId(), unphased_ptr->populationId());

    if (not InbreedingAnalysis::populationInbreeding(genome_ptr,
                                                     unphased_ptr,
                                                     *diploid_ptr,
                                                     *ped_data_ptr,
                                                     parameters)) {

      ExecEnv::log().error("InbreedingAnalysis::Inbreeding; problem with population inbreeding analysis");
      return false;

    }

  } else if (unphased_ptr) {

    if (not InbreedingAnalysis::syntheticInbreeding(unphased_ptr, parameters)) {

      ExecEnv::log().error("InbreedingAnalysis::Inbreeding;problem with synthetic inbreeding analysis");
      return false;

    }

  } else {

    ExecEnv::log().error("InbreedingAnalysis::Inbreeding, Failed to create Filtered Population and Joined Populations");
    return false;

  }

  return true;

}



// Calculate the Population Inbreeding Coefficient
bool kgl::InbreedingAnalysis::populationInbreeding(std::shared_ptr<const GenomeReference> genome_ptr,
                                                   std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const DiploidPopulation& diploid_population,
                                                   const GenomePEDData& ped_data,
                                                   InbreedingParameters& parameters) {

  static const ContigOffset_t sampling_distance = 300;
  static const size_t locii_count = 10000;

  // Filter out any variants that did not pass VCF filters (otherwise we get duplicate variants).
  unphased_ptr = unphased_ptr->filterVariants(PassFilter());

  // check that unphased population onlu has 1 genome
  if (unphased_ptr->getMap().size() != 1) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding; Unphased Population: {} has unexpected Genome count: {}",
                         unphased_ptr->populationId(), unphased_ptr->getMap().size());
    return false;

  }

  // check that unphased only has 1 contig
  auto [genone_id, contig_map] = *unphased_ptr->getMap().begin();
  if (contig_map->getMap().size() != 1) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding; Unphased Population: {} Genome: {} has more than 1 contig: {}",
                         unphased_ptr->populationId(), genone_id, contig_map->getMap().size());
    return false;

  }

  // Get the size of the contig.
  auto [contig_id, contig_ptr] = *contig_map->getMap().begin();
  auto contig_opt = genome_ptr->getContigSequence(contig_id);

  if (not contig_opt) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding; Unphased Population: {} Genome: {} not found in the reference genome",
                         unphased_ptr->populationId(), genone_id);
    return false;
  }

  ContigOffset_t contig_size = contig_opt.value()->contigSize();
  ContigOffset_t lower_window = 0;

  parameters.lociiArguments().lociiCount(locii_count);
  parameters.lociiArguments().lowerOffset(lower_window);
  parameters.lociiArguments().lociiSpacing(sampling_distance);

  std::vector<ContigOffset_t> locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                                                VariantDatabaseRead::SUPER_POP_ALL_,
                                                                                parameters.lociiArguments());
  ContigOffset_t upper_window = locii_vector.back();
  parameters.lociiArguments().upperOffset(upper_window);

  while (lower_window < (contig_size - 1000) and locii_vector.size() >= 1000) {

    std::stringstream ss;

    ss << "_" << lower_window << "_" << upper_window << "_";
    std::string window_file = parameters.outputFile() + ss.str();


    populationInbreedingSample(unphased_ptr,
                               diploid_population,
                               ped_data,
                               window_file,
                               parameters);

    lower_window = upper_window;
    parameters.lociiArguments().lowerOffset(lower_window);
    locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                      VariantDatabaseRead::SUPER_POP_ALL_,
                                                      parameters.lociiArguments());
    upper_window = locii_vector.back();
    parameters.lociiArguments().upperOffset(upper_window);

  }

  return true;

}


bool kgl::InbreedingAnalysis::populationInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                          const DiploidPopulation& diploid_population,
                                                          const GenomePEDData& ped_data,
                                                          const std::string& output_file_name,
                                                          const InbreedingParameters& parameters) {

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = output_file_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreedingSample; could not open output file: {}", file_name_ext);
    return false;

  }

  outfile << diploid_population.populationId()
          << DELIMITER_ << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing()
          << DELIMITER_ << "LowerOffset:" << parameters.lociiArguments().lowerOffset()
          << DELIMITER_ << "UpperOffset:" << parameters.lociiArguments().upperOffset() << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr, parameters.lociiArguments());

  processResults(contig_locus_map, diploid_population, ped_data, outfile, parameters);

  return true;

}



bool kgl::InbreedingAnalysis::processResults( const ContigLocusMap& contig_locus_map,
                                              const DiploidPopulation& diploid_population,
                                              const GenomePEDData& ped_data,
                                              std::ostream& outfile,
                                              const InbreedingParameters& parameters) {

  // Retrieve the algorithm function object.
  auto algorithm_opt = InbreedingCalculation::namedAlgorithm(parameters.inbreedingAlgorthim());

  if (not algorithm_opt) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Inbreeding algorithm not found: {}", parameters.inbreedingAlgorthim());
    return false;

  }

  for (auto const& [genome_contig_id, locus_map] : contig_locus_map) {

    // Use a thread pool to calculate inbreeding and relatedness.
    ThreadPool thread_pool;
    std::vector<std::future<LocusResults>> future_vector;
    for (auto const&[genome_id, genome_ptr] : diploid_population.getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto result = ped_data.getMap().find(genome_id);
        if (result == ped_data.getMap().end()) {

          ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Genome sample: {} does not have a PED record", genome_id);
          continue;

        }
        auto const& [sample_id, ped_record] = *result;
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<LocusResults> future = thread_pool.enqueueTask(algorithm_opt.value(),
                                                                   genome_id,
                                                                   contig_opt.value(),
                                                                   super_pop_id,
                                                                   locus_list,
                                                                   parameters);
        future_vector.push_back(std::move(future));

      }

    }

    // Retrieve the thread results into a map.
    ResultsMap genome_results_map;
    for (auto &future : future_vector) {

      auto locus_results = future.get();
      ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
      genome_results_map[locus_results.genome] = locus_results;


    }

    if (not genome_results_map.empty()) {

      writeResults(genome_contig_id, genome_results_map, ped_data, outfile);

    }

  }

  return true;

}


bool kgl::InbreedingAnalysis::writeResults( const ContigId_t& contig_id,
                                            const ResultsMap& genome_results_map,
                                            const GenomePEDData& ped_data,
                                            std::ostream& outfile) {


  outfile << contig_id
          << DELIMITER_ << "Population"
          << DELIMITER_ << "Description"
          << DELIMITER_ << "SuperPopulation"
          << DELIMITER_ << "Description"
          << DELIMITER_ << "MajorHet"
          << DELIMITER_ << "MajorHetFreq"
          << DELIMITER_ << "MinorHet"
          << DELIMITER_ << "MinorHetFreq"
          << DELIMITER_ << "MajorHom"
          << DELIMITER_ << "MajorHomFreq"
          << DELIMITER_ << "MinorHom"
          << DELIMITER_ << "MinorHomFreq"
          << DELIMITER_ << "Het/Hom"
          << DELIMITER_ << "TotalLoci"
          << DELIMITER_ << "CalcInbreed" << '\n';

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto result = ped_data.getMap().find(genome_id);

    if (result == ped_data.getMap().end()) {

      ExecEnv::log().error("InbreedingAnalysis::writeResults, Genome sample: {} does not have a PED record", genome_id);
      continue;

    }

    auto const& [sample_id, ped_record] = *result;

    double het_hom_ratio;
    if (locus_results.minor_homo_count > 0) {

      size_t total_het = locus_results.major_hetero_count + locus_results.minor_hetero_count;
      het_hom_ratio = static_cast<double>(total_het) / static_cast<double>(locus_results.minor_homo_count);

    } else {

      het_hom_ratio = 0.0;

    }

    outfile << genome_id << DELIMITER_;
    outfile << ped_record.population() << DELIMITER_;
    outfile << ped_record.populationDescription() << DELIMITER_;
    outfile << ped_record.superPopulation() << DELIMITER_;
    outfile << ped_record.superDescription() << DELIMITER_;
    outfile << locus_results.major_hetero_count << DELIMITER_;
    outfile << locus_results.major_hetero_freq << DELIMITER_;
    outfile << locus_results.minor_hetero_count << DELIMITER_;
    outfile << locus_results.minor_hetero_freq << DELIMITER_;
    outfile << locus_results.major_homo_count << DELIMITER_;
    outfile << locus_results.major_homo_freq << DELIMITER_;
    outfile << locus_results.minor_homo_count << DELIMITER_;
    outfile << locus_results.minor_homo_freq << DELIMITER_;
    outfile << het_hom_ratio << DELIMITER_;
    outfile << locus_results.total_allele_count << DELIMITER_;
    outfile << locus_results.inbred_allele_sum << DELIMITER_;
    outfile << '\n';

  }

  outfile.flush();

  return true;

}


// Calculate the Synthetic Population Inbreeding Coefficient
bool kgl::InbreedingAnalysis::syntheticInbreeding(  std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                    InbreedingParameters& parameters) {


  static const ContigOffset_t sampling_distance = 10000;
  // Filter out any variants that did not pass VCF filters (otherwise we get duplicate variants).
  unphased_ptr = unphased_ptr->filterVariants(PassFilter());

  parameters.lociiArguments().lociiSpacing(sampling_distance);

  syntheticInbreedingSample( unphased_ptr, parameters);


  return true;

}


bool kgl::InbreedingAnalysis::syntheticInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                         InbreedingParameters& parameters) {

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = parameters.outputFile() + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::syntheticInbreedingSample; could not open output file: {}", file_name_ext);
    return false;

  }

  outfile << "Synthetic"
          << DELIMITER_ << "Min_AF:" << parameters.lociiArguments().minAlleleFrequency()
          << DELIMITER_ << "Max_AF:" << parameters.lociiArguments().maxAlleleFrequency()
          << DELIMITER_ << "Spacing:" << parameters.lociiArguments().lociiSpacing() << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr, parameters.lociiArguments());  // spacing

  processSynResults(contig_locus_map, outfile, parameters);

  return true;

}



bool kgl::InbreedingAnalysis::processSynResults( const ContigLocusMap& contig_locus_map,
                                                 std::ostream& outfile,
                                                 InbreedingParameters& parameters) {

  // Retrieve the algorithm function object.
  auto algorithm_opt = InbreedingCalculation::namedAlgorithm(parameters.inbreedingAlgorthim());

  if (not algorithm_opt) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Inbreeding algorithm not found: {}", parameters.inbreedingAlgorthim());
    return false;

  }

  // For each contig.
  for (auto const& [genome_contig_id, locus_map] : contig_locus_map) {

    // Use a thread pool to calculate inbreeding and relatedness.
    for (auto const&[super_pop_id, locus_list] : locus_map) {

      std::shared_ptr<const DiploidPopulation> population = InbreedSynthetic::generateSyntheticPopulation( MIN_INBREEDING_COEFICIENT,
                                                                                                          MAX_INBREEDING_COEFICIENT,
                                                                                                          STEP_INBREEDING_COEFICIENT,
                                                                                                          super_pop_id,
                                                                                                          *locus_list,
                                                                                                          parameters.lociiArguments());

      ThreadPool thread_pool;
      std::vector<std::future<LocusResults>> future_vector;
      for (auto const&[genome_id, genome_ptr] : population->getMap()) {

        auto contig_opt = genome_ptr->getContig(genome_contig_id);

        if (contig_opt) {


          std::future<LocusResults> future = thread_pool.enqueueTask(algorithm_opt.value(),
                                                                     genome_id,
                                                                     contig_opt.value(),
                                                                     super_pop_id,
                                                                     locus_list,
                                                                     parameters);
          future_vector.push_back(std::move(future));

        }

      }

      // Retrieve the thread results into a map.
      ResultsMap genome_results_map;
      for (auto &future : future_vector) {

        auto locus_results = future.get();
        ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
        genome_results_map[locus_results.genome] = locus_results;


      }

      if (not genome_results_map.empty()) {

        writeSynResults(genome_contig_id, genome_results_map, outfile);

      }

    } // for locus

  } // for contig.

  return true;

}


bool kgl::InbreedingAnalysis::writeSynResults(const ContigId_t& contig_id,
                                              const std::map<GenomeId_t, LocusResults>& genome_results_map,
                                              std::ostream& outfile) {

  outfile << contig_id << DELIMITER_
          << "SuperPop" << DELIMITER_
          << "Inbreeding" << DELIMITER_
          << "MajorHet" << DELIMITER_
          << "MajorHetFreq" << DELIMITER_
          << "MinorHet" << DELIMITER_
          << "MinorHetFreq" << DELIMITER_
          << "MajorHom" << DELIMITER_
          << "MajorHomFreq" << DELIMITER_
          << "MinorHom" << DELIMITER_
          << "MinorHomFreq" << DELIMITER_
          << "Het/Hom"<< DELIMITER_
          << "TotalLoci" << DELIMITER_
          << "CalcInbreed" << '\n';

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto const [valid_value, inbreeding] = InbreedSynthetic::generateInbreeding(genome_id);

    double het_hom_ratio;
    if (locus_results.minor_homo_count > 0) {

      size_t total_het = locus_results.major_hetero_count + locus_results.minor_hetero_count;
      het_hom_ratio = static_cast<double>(total_het) / static_cast<double>(locus_results.minor_homo_count);

    } else {

      het_hom_ratio = 0.0;

    }

    outfile << genome_id << DELIMITER_;
    outfile << genome_id.substr(0, genome_id.find_first_of("_")) << DELIMITER_;

    if (valid_value) {

      outfile << inbreeding << DELIMITER_;

    } else {

      outfile << "Invalid" << DELIMITER_;

    }
    outfile << locus_results.major_hetero_count << DELIMITER_;
    outfile << locus_results.major_hetero_freq << DELIMITER_;
    outfile << locus_results.minor_hetero_count << DELIMITER_;
    outfile << locus_results.minor_hetero_freq << DELIMITER_;
    outfile << locus_results.major_homo_count << DELIMITER_;
    outfile << locus_results.major_homo_freq << DELIMITER_;
    outfile << locus_results.minor_homo_count << DELIMITER_;
    outfile << locus_results.minor_homo_freq << DELIMITER_;
    outfile << het_hom_ratio << DELIMITER_;
    outfile << locus_results.total_allele_count << DELIMITER_;
    outfile << locus_results.inbred_allele_sum << DELIMITER_;
    outfile << '\n';

  }

  outfile.flush();

  return true;

}

