//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_inbreed_aux.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Calculate the Population Inbreeding Coefficient
bool kgl::InbreedingAnalysis::populationInbreeding(std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const DiploidPopulation& diploid_population,
                                                   const GenomePEDData& ped_data,
                                                   const std::string& output_file_name) {


  populationInbreedingSample( unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_all",
                              0.0,
                              1.0,
                              10);

  populationInbreedingSample( unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_100_01",
                              0.01,
                              1.0,
                              10);

  populationInbreedingSample( unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_05_01",
                              0.01,
                              0.05,
                              10);

  populationInbreedingSample( unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_20_05",
                              0.05,
                              0.2,
                              10);

  populationInbreedingSample( unphased_ptr,
                              diploid_population,
                              ped_data,
                              output_file_name,
                              "pop_100_20",
                              0.2,
                              1.0,
                              10);


  return true;

}


bool kgl::InbreedingAnalysis::populationInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                          const DiploidPopulation& diploid_population,
                                                          const GenomePEDData& ped_data,
                                                          const std::string& output_file_name,
                                                          const std::string& sample_name,
                                                          double min,
                                                          double max,
                                                          ContigOffset_t  spacing) {

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = output_file_name + sample_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile << diploid_population.populationId() << DELIMITER_ << "Min_AF:" << DELIMITER_ << min
           << DELIMITER_ << "Max_AF:" << max << DELIMITER_ << "Spacing:" << spacing << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                                           min,   // min af
                                                                           max,   // max af
                                                                           spacing);  // spacing

  processResults(contig_locus_map, diploid_population, ped_data, outfile);

  return true;

}



bool kgl::InbreedingAnalysis::processResults( const ContigLocusMap& contig_locus_map,
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

        std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingCalculation::processLogLikelihood,
                                                                   genome_id,
                                                                   contig_opt.value(),
                                                                   InbreedSampling::lookupSuperPopulationField(super_pop_id),
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
          << "HetCount" << DELIMITER_ << "HomCount" << DELIMITER_ << "Het/Hom"<< DELIMITER_
          << "TotalLoci" << DELIMITER_ << "CalcInbreed\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto result = ped_data.getMap().find(genome_id);

    if (result == ped_data.getMap().end()) {

      ExecEnv::log().error("InbreedingAnalysis::writeResults, Genome sample: {} does not have a PED record", genome_id);
      continue;

    }

    auto const& [sample_id, ped_record] = *result;

    double het_hom_ratio = (locus_results.homo_count > 0 ? static_cast<double>(locus_results.hetero_count)/static_cast<double>(locus_results.homo_count) : 0.0);

    outfile << genome_id << DELIMITER_;
    outfile << ped_record.population() << DELIMITER_;
    outfile << ped_record.populationDescription() << DELIMITER_;
    outfile << ped_record.superPopulation() << DELIMITER_;
    outfile << ped_record.superDescription() << DELIMITER_;
    outfile << locus_results.hetero_count << DELIMITER_;
    outfile << locus_results.homo_count << DELIMITER_;
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
                                                    const std::string& output_file_name) {

/*
  syntheticInbreedingSample( unphased_ptr,
                             output_file_name,
                             "syn_all",
                             0.0,
                             1.0,
                             10);

  syntheticInbreedingSample( unphased_ptr,
                             output_file_name,
                             "syn_100_01",
                             0.01,
                             1.0,
                             10);

  syntheticInbreedingSample( unphased_ptr,
                             output_file_name,
                             "syn_05_01",
                             0.01,
                             0.05,
                             10);

  syntheticInbreedingSample( unphased_ptr,
                             output_file_name,
                             "syn_20_05",
                             0.05,
                             0.2,
                             10);

  syntheticInbreedingSample( unphased_ptr,
                             output_file_name,
                             "syn_100_20",
                             0.2,
                             1.0,
                             10);
*/
  return true;

}


bool kgl::InbreedingAnalysis::syntheticInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                          const std::string& output_file_name,
                                                          const std::string& sample_name,
                                                          double min,
                                                          double max,
                                                          ContigOffset_t  spacing) {

  // Open the output file.
  std::ofstream outfile;
  std::string file_name_ext = output_file_name + sample_name + FILE_EXT_;
  outfile.open(file_name_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile << "Synthetic" << DELIMITER_ << "Min_AF:" << DELIMITER_ << min
          << DELIMITER_ << "Max_AF:" << max << DELIMITER_ << "Spacing:" << spacing << '\n';

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr,
                                                                           min,   // min af
                                                                           max,   // max af
                                                                           spacing);  // spacing

  processSynResults(contig_locus_map, outfile);

  return true;

}



bool kgl::InbreedingAnalysis::processSynResults( const ContigLocusMap& contig_locus_map,
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


          std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingCalculation::processLogLikelihood,
                                                                     genome_id,
                                                                     contig_opt.value(),
                                                                     InbreedSampling::lookupSuperPopulationField(super_pop_id),
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


  outfile << contig_id << DELIMITER_ << "SuperPop" << DELIMITER_ << "Inbreeding"
          << DELIMITER_ << "HetCount" << DELIMITER_ << "HomCount"
          << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "CalcInbreed\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {


    auto const [valid_value, inbreeding] = InbreedSampling::generateInbreeding(genome_id);
    double het_hom_ratio = (locus_results.homo_count > 0 ? static_cast<double>(locus_results.hetero_count)/static_cast<double>(locus_results.homo_count) : 0.0);

    outfile << genome_id << DELIMITER_;
    outfile << genome_id.substr(0, genome_id.find_first_of("_")) << DELIMITER_;

    if (valid_value) {

      outfile << inbreeding << DELIMITER_;

    } else {

      outfile << "Invalid" << DELIMITER_;

    }
    outfile << locus_results.hetero_count << DELIMITER_;
    outfile << locus_results.homo_count << DELIMITER_;
    outfile << het_hom_ratio << DELIMITER_;
    outfile << locus_results.total_allele_count << DELIMITER_;
    outfile << locus_results.inbred_allele_sum << DELIMITER_;
    outfile << '\n';

  }

  outfile.flush();

  return true;

}

