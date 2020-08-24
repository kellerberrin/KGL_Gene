//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_inbreed_aux.h"


#include <fstream>

namespace kgl = kellerberrin::genome;


// Calculate the HetHom Ratio
bool kgl::InbreedingAnalysis::populationInbreeding(const GenomeReference& genome_GRCh38,
                                                   const UnphasedPopulation& unphased_population,
                                                   const DiploidPopulation& diploid_population,
                                                   const GenomePEDData& ped_data,
                                                   const std::string& output_file_name) {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38.getMap()) {

    // Generate the super population allele locus lists
    LocusMap locus_map = InbreedSampling::getPopulationLocus(unphased_population, genome_contig_id);

    // Use a thread pool to calculate inbreeding and relatedness.
    ThreadPool thread_pool;
    std::vector<std::future<LocusResults>> future_vector;
    for (auto const&[genome_id, genome_ptr] : diploid_population.getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto result = ped_data.getMap().find(genome_id);
        if (result == ped_data.getMap().end()) {

          ExecEnv::log().error("InbreedingAnalysis::hetHomRatioLocu, Genome sample: {} does not have a PED record", genome_id);
          continue;

        }
        auto const& [sample_id, ped_record] = *result;
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("InbreedingAnalysis::hetHomRatioLocu, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingAnalysis::processRitlandLocus,
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

      writeResults(genome_contig_id, genome_results_map, ped_data, output_file_name);

    }

  }

  return true;

}



bool kgl::InbreedingAnalysis::writeResults( const ContigId_t& contig_id,
                                            const ResultsMap& genome_results_map,
                                            const GenomePEDData& ped_data,
                                            const std::string& output_file_name) {

  // Append the results.
  std::ofstream outfile;
  std::string file_ext = output_file_name + FILE_EXT_;
  outfile.open(file_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile << contig_id << DELIMITER_ << "Population" << DELIMITER_
          << "Description" << DELIMITER_ << "SuperPopulation" << DELIMITER_ << "Description" << DELIMITER_
          << "HetCount" << DELIMITER_ << "HomCount" << DELIMITER_ << "Het/Hom"<< DELIMITER_
          << "TotalLoci" << DELIMITER_ << "Ritland\n";

  for (auto const& [genome_id, locus_results] : genome_results_map) {

    auto result = ped_data.getMap().find(genome_id);

    if (result == ped_data.getMap().end()) {

      ExecEnv::log().error("InbreedingAnalysis::hetHomRatio, Genome sample: {} does not have a PED record", genome_id);
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


// Calculate the HetHom Ratio
bool kgl::InbreedingAnalysis::syntheticInbreeding(  const GenomeReference& genome_GRCh38,
                                                    const UnphasedPopulation& unphased_population,
                                                    const std::string& output_file_name) {

  for (auto const& [genome_contig_id, genome_contig_ptr] : genome_GRCh38.getMap()) {

    // Generate the super population allele locus lists
    LocusMap locus_map = InbreedSampling::getPopulationLocus(unphased_population, genome_contig_id);

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


          std::future<LocusResults> future = thread_pool.enqueueTask(&InbreedingAnalysis::processRitlandLocus,
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

        syntheticResults(genome_contig_id, genome_results_map, output_file_name);

      }

    } // for locus

  } // for contig.

  return true;

}



bool kgl::InbreedingAnalysis::syntheticResults( const ContigId_t& contig_id,
                                                const std::map<GenomeId_t, LocusResults>& genome_results_map,
                                                const std::string& output_file_name) {

  // Append the results.
  std::ofstream outfile;
  std::string file_ext = output_file_name + FILE_EXT_;
  outfile.open(file_ext, std::ofstream::out |  std::ofstream::trunc);

  outfile << contig_id << DELIMITER_ << "SuperPop" << DELIMITER_ << "Inbreeding"
          << DELIMITER_ << "HetCount" << DELIMITER_ << "HomCount"
          << DELIMITER_ << "Het/Hom"<< DELIMITER_ << "TotalLoci" << DELIMITER_ << "Ritland\n";

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

