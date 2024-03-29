//
// Created by kellerberrin on 13/8/20.
//


#include "kel_workflow_threads.h"
#include "kga_analysis_inbreed_diploid.h"
#include "kga_analysis_inbreed_locus.h"
#include "kga_analysis_inbreed_output.h"

#include <iostream>


namespace kga = kellerberrin::genome::analysis;


// Calculate the Population Inbreeding Coefficient
bool kga::InbreedingAnalysis::populationInbreeding(std::shared_ptr<const PopulationDB> unphased_ptr,
                                                   const PopulationDB& diploid_population,
                                                   const HsGenomeGenealogyData& ped_data,
                                                   InbreedParamOutput& param_output) {



  // check that unphased population onlu has 1 genome
  if (unphased_ptr->getMap().size() != 1) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding; Unphased Population: {} has unexpected Genome count: {}",
                         unphased_ptr->populationId(), unphased_ptr->getMap().size());
    return false;

  }

  // check that unphased only has 1 contig_ref_ptr
  auto [genone_id, contig_map] = *unphased_ptr->getMap().begin();
  if (contig_map->getMap().size() != 1) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding; Unphased Population: {} Genome: {} has more than 1 contig_ref_ptr: {}",
                         unphased_ptr->populationId(), genone_id, contig_map->getMap().size());
    return false;

  }

  // Get the size of the contig_ref_ptr.
  auto [contig_id, contig_ptr] = *contig_map->getMap().begin();

  InbreedingParameters local_params = param_output.getParameters();
  std::vector<ContigOffset_t> locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                                                FrequencyDatabaseRead::SUPER_POP_ALL_,
                                                                                local_params.lociiArguments());
  local_params.lociiArguments().upperOffset(locii_vector.back());

  while ( local_params.lociiArguments().upperOffset() < param_output.getParameters().lociiArguments().upperOffset()
          and locii_vector.size() >= 100) {

    ResultsMap results_map = populationInbreedingSample( unphased_ptr,
                                                         diploid_population,
                                                         ped_data,
                                                         local_params);

    // Generate a string to identify these results.
    std::string result_ident = InbreedingResultColumn::generateIdent( contig_id,
                                                                      local_params.lociiArguments().lowerOffset(),
                                                                      local_params.lociiArguments().upperOffset());

    // Store the inbreeding results.
    param_output.addColumn(InbreedingResultColumn(result_ident, results_map));

    local_params.lociiArguments().lowerOffset(local_params.lociiArguments().upperOffset());
    locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                      FrequencyDatabaseRead::SUPER_POP_ALL_,
                                                      local_params.lociiArguments());
    local_params.lociiArguments().upperOffset(locii_vector.back());

  }

  return true;

}



kga::ResultsMap kga::InbreedingAnalysis::populationInbreedingSample( std::shared_ptr<const PopulationDB> unphased_ptr,
                                                                     const PopulationDB& diploid_population,
                                                                     const HsGenomeGenealogyData& ped_data,
                                                                     const InbreedingParameters& parameters) {

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr, parameters.lociiArguments());

  ResultsMap results_map = processResults(contig_locus_map, diploid_population, ped_data, parameters);

  return results_map;

}



kga::ResultsMap kga::InbreedingAnalysis::processResults( const ContigLocusMap& contig_locus_map,
                                                         const PopulationDB& diploid_population,
                                                         const HsGenomeGenealogyData& ped_data,
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
    WorkflowThreads thread_pool(WorkflowThreads::defaultThreads());
    std::vector<std::future<LocusResults>> future_vector;
    for (auto const&[genome_id, genome_ptr] : diploid_population.getMap()) {

      auto contig_opt = genome_ptr->getContig(genome_contig_id);

      if (contig_opt) {

        auto record_opt = ped_data.getGenomeGenealogyRecord(genome_id);
        if (not record_opt) {

          ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Genome sample: {} does not have a PED record", genome_id);
          continue;

        }
        auto const ped_record = record_opt.value();
        auto locus_result = locus_map.find(ped_record.superPopulation());
        if (locus_result == locus_map.end()) {

          ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Locus set not found for super population: {}", ped_record.superPopulation());
          continue;

        }
        auto const& [super_pop_id, locus_list] = *locus_result;

        std::future<LocusResults> future = thread_pool.enqueueFuture(algorithm_opt.value(),
                                                                     genome_id,
                                                                     contig_opt.value(),
                                                                     super_pop_id,
                                                                     locus_list);
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


