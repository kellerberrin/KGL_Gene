//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_inbreed_aux.h"
#include "kgl_filter.h"


#include <fstream>

namespace kgl = kellerberrin::genome;



bool kgl::InbreedingAnalysis::Inbreeding( InbreedingAlgorithm algorithm,
                                          std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                          std::shared_ptr<const DiploidPopulation> diploid_ptr,
                                          std::shared_ptr<const GenomePEDData> ped_data_ptr,
                                          const std::string& output_file_name) {


  if (diploid_ptr and unphased_ptr and ped_data_ptr) {

    ExecEnv::log().info("Filtered Population : {}  and Joined Population: {}  both active",
                        diploid_ptr->populationId(), unphased_ptr->populationId());

    if (not InbreedingAnalysis::populationInbreeding(algorithm,
                                                     unphased_ptr,
                                                     *diploid_ptr,
                                                     *ped_data_ptr,
                                                     output_file_name)) {

      ExecEnv::log().error("InbreedingAnalysis::Inbreeding; problem with population inbreeding analysis");
      return false;

    }

    if (not InbreedingAnalysis::syntheticInbreeding(algorithm,
                                                    unphased_ptr,
                                                    output_file_name)) {

      ExecEnv::log().error("InbreedingAnalysis::Inbreeding; problem with synthetic inbreeding analysis");
      return false;

    }

  } else if (unphased_ptr){

    if (not InbreedingAnalysis::syntheticInbreeding(algorithm,
                                                    unphased_ptr,
                                                    output_file_name)) {

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
bool kgl::InbreedingAnalysis::populationInbreeding(InbreedingAlgorithm algorithm,
                                                   std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const DiploidPopulation& diploid_population,
                                                   const GenomePEDData& ped_data,
                                                   const std::string& output_file_name) {

  const ContigOffset_t sampling_distance = 10000;

  // Filter out any variants that did not pass VCF filters (otherwise we get duplicate variants).
  unphased_ptr = unphased_ptr->filterVariants(PassFilter());

  populationInbreedingSample( algorithm,
                              unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_all",
                              0.0,
                              1.0,
                              sampling_distance);

  populationInbreedingSample( algorithm,
                              unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_100_01",
                              0.01,
                              1.0,
                              sampling_distance);

  populationInbreedingSample( algorithm,
                              unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_05_01",
                              0.01,
                              0.05,
                              sampling_distance);

  populationInbreedingSample( algorithm,
                              unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_20_05",
                              0.05,
                              0.2,
                              sampling_distance);

  populationInbreedingSample( algorithm,
                              unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_100_20",
                              0.2,
                              1.0,
                              sampling_distance);


  return true;

}


bool kgl::InbreedingAnalysis::populationInbreedingSample( InbreedingAlgorithm algorithm,
                                                          std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                          const DiploidPopulation& diploid_population,
                                                          const GenomePEDData& ped_data,
                                                          const std::string& output_file_name,
                                                          const std::string& sample_name,
                                                          double allele_frequency_min,
                                                          double allele_frequency_max,
                                                          ContigOffset_t  spacing) {

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = output_file_name + sample_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::syntheticInbreedingSample; could not open output file: {}", file_name_ext);
    return false;

  }

  outfile << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << allele_frequency_min
          << DELIMITER_ << "Max_AF:" << allele_frequency_max << DELIMITER_ << "Spacing:" << spacing << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                                           allele_frequency_min,   // min af
                                                                           allele_frequency_max,   // max af
                                                                           spacing);  // spacing

  processResults(algorithm, contig_locus_map, diploid_population, ped_data, outfile);

  return true;

}



bool kgl::InbreedingAnalysis::processResults( InbreedingAlgorithm algorithm,
                                              const ContigLocusMap& contig_locus_map,
                                              const DiploidPopulation& diploid_population,
                                              const GenomePEDData& ped_data,
                                              std::ostream& outfile) {

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

        std::future<LocusResults> future = thread_pool.enqueueTask(algorithm,
                                                                   genome_id,
                                                                   contig_opt.value(),
                                                                   super_pop_id,
                                                                   locus_list);
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


  outfile << contig_id << DELIMITER_ << "Population" << DELIMITER_
          << "Description" << DELIMITER_ << "SuperPopulation" << DELIMITER_ << "Description" << DELIMITER_
          << "MajorHet" << DELIMITER_ << "MajorHetFreq" << DELIMITER_ << "MinorHet" << DELIMITER_ << "MinorHetFreq"
          << DELIMITER_ << "MajorHom" << DELIMITER_ << "MajorHomFreq" << DELIMITER_ << "MinorHom" << DELIMITER_
          << "MinorHomFreq" << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "CalcInbreed\n";

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
bool kgl::InbreedingAnalysis::syntheticInbreeding(  InbreedingAlgorithm algorithm,
                                                    std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                    const std::string& output_file_name) {

  const ContigOffset_t sampling_distance = 10000;

  // Filter out any variants that did not pass VCF filters (otherwise we get duplicate variants).
  unphased_ptr = unphased_ptr->filterVariants(PassFilter());

  syntheticInbreedingSample( algorithm,
                             unphased_ptr,
                             output_file_name,
                             "syn_all",
                             0.0,
                             1.0,
                             sampling_distance);

  syntheticInbreedingSample( algorithm,
                             unphased_ptr,
                             output_file_name,
                             "syn_100_01",
                             0.01,
                             1.0,
                             sampling_distance);

  syntheticInbreedingSample( algorithm,
                             unphased_ptr,
                             output_file_name,
                             "syn_05_01",
                             0.01,
                             0.05,
                             sampling_distance);

  syntheticInbreedingSample( algorithm,
                             unphased_ptr,
                             output_file_name,
                             "syn_20_05",
                             0.05,
                             0.2,
                             sampling_distance);

  syntheticInbreedingSample( algorithm,
                             unphased_ptr,
                             output_file_name,
                             "syn_100_20",
                             0.2,
                             1.0,
                             sampling_distance);

  return true;

}


bool kgl::InbreedingAnalysis::syntheticInbreedingSample( InbreedingAlgorithm algorithm,
                                                         std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                         const std::string& output_file_name,
                                                         const std::string& sample_name,
                                                         double allele_frequency_min,
                                                         double allele_frequency_max,
                                                         ContigOffset_t  spacing) {

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = output_file_name + sample_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  if (not outfile.good()) {

    ExecEnv::log().error("InbreedingAnalysis::syntheticInbreedingSample; could not open output file: {}", file_name_ext);
    return false;

  }

  outfile << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << allele_frequency_min
          << DELIMITER_ << "Max_AF:" << allele_frequency_max << DELIMITER_ << "Spacing:" << spacing << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                                           allele_frequency_min,   // min af
                                                                           allele_frequency_max,   // max af
                                                                           spacing);  // spacing

  processSynResults(algorithm, contig_locus_map, outfile);

  return true;

}



bool kgl::InbreedingAnalysis::processSynResults( InbreedingAlgorithm algorithm,
                                                 const ContigLocusMap& contig_locus_map,
                                                 std::ostream& outfile) {

  // For each contig.
  for (auto const& [genome_contig_id, locus_map] : contig_locus_map) {

    // Use a thread pool to calculate inbreeding and relatedness.
    for (auto const&[super_pop_id, locus_list] : locus_map) {

      std::shared_ptr<const DiploidPopulation> population = InbreedSampling::generateSyntheticPopulation( MIN_INBREEDING_COEFICIENT,
                                                                                                          MAX_INBREEDING_COEFICIENT,
                                                                                                          STEP_INBREEDING_COEFICIENT,
                                                                                                          super_pop_id,
                                                                                                          *locus_list);

      ThreadPool thread_pool;
      std::vector<std::future<LocusResults>> future_vector;
      for (auto const&[genome_id, genome_ptr] : population->getMap()) {

        auto contig_opt = genome_ptr->getContig(genome_contig_id);

        if (contig_opt) {


          std::future<LocusResults> future = thread_pool.enqueueTask(algorithm,
                                                                     genome_id,
                                                                     contig_opt.value(),
                                                                     super_pop_id,
                                                                     locus_list);
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


  outfile << contig_id << DELIMITER_ << "SuperPop" << DELIMITER_ << "Inbreeding" << DELIMITER_
          << "MajorHet" << DELIMITER_ << "MajorHetFreq" << DELIMITER_ << "MinorHet" << DELIMITER_ << "MinorHetFreq"
          << DELIMITER_ << "MajorHom" << DELIMITER_ << "MajorHomFreq" << DELIMITER_ << "MinorHom" << DELIMITER_
          << "MinorHomFreq" << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "CalcInbreed\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {


    auto const [valid_value, inbreeding] = InbreedSampling::generateInbreeding(genome_id);

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

