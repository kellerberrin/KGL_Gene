//
// Created by kellerberrin on 13/8/20.
//


#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed.h"
#include "kgl_analysis_mutation_inbreed_locus.h"
#include "kgl_analysis_mutation_inbreed_output.h"

#include <iostream>


namespace kgl = kellerberrin::genome;




// Calculate the Population Inbreeding Coefficient
bool kgl::InbreedingAnalysis::populationInbreeding(std::shared_ptr<const PopulationVariant> unphased_ptr,
                                                   const PopulationVariant& diploid_population,
                                                   const GenomePEDData& ped_data,
                                                   const InbreedingParameters& parameters,
                                                   InbreedingOutputResults& results) {



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

  InbreedingParameters local_params = parameters;
  std::vector<ContigOffset_t> locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                                                FrequencyDatabaseRead::SUPER_POP_ALL_,
                                                                                local_params.lociiArguments());
  local_params.lociiArguments().upperOffset(locii_vector.back());

  while (local_params.lociiArguments().upperOffset() < parameters.lociiArguments().upperOffset() and locii_vector.size() >= 100) {

    ResultsMap results_map = populationInbreedingSample( unphased_ptr,
                                                         diploid_population,
                                                         ped_data,
                                                         local_params);

    results.insertResults({local_params, results_map});

    local_params.lociiArguments().lowerOffset(local_params.lociiArguments().upperOffset());
    locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                      FrequencyDatabaseRead::SUPER_POP_ALL_,
                                                      local_params.lociiArguments());
    local_params.lociiArguments().upperOffset(locii_vector.back());

  }

  return true;

}



kgl::ResultsMap kgl::InbreedingAnalysis::populationInbreedingSample( std::shared_ptr<const PopulationVariant> unphased_ptr,
                                                                     const PopulationVariant& diploid_population,
                                                                     const GenomePEDData& ped_data,
                                                                     const InbreedingParameters& parameters) {

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr, parameters.lociiArguments());

  ResultsMap results_map = processResults(contig_locus_map, diploid_population, ped_data, parameters);

  return results_map;

}



kgl::ResultsMap kgl::InbreedingAnalysis::processResults( const ContigLocusMap& contig_locus_map,
                                                         const PopulationVariant& diploid_population,
                                                         const GenomePEDData& ped_data,
                                                         const InbreedingParameters& parameters) {

  ResultsMap results_map;
  // Retrieve the algorithm function object.
  auto algorithm_opt = InbreedingCalculation::namedAlgorithm(parameters.inbreedingAlgorthim());

  if (not algorithm_opt) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Inbreeding algorithm not found: {}", parameters.inbreedingAlgorthim());
    return results_map;

  }

  for (auto const& [genome_contig_id, locus_map] : contig_locus_map) {

    // Use a thread pool to calculate inbreeding and relatedness.
    ThreadPool thread_pool(ThreadPool::hardwareThreads());
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
    for (auto &future : future_vector) {

      auto locus_results = future.get();
      ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
      results_map[locus_results.genome] = locus_results;

    }

  } // all contigs.

  return results_map;

}


