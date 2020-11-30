//
// Created by kellerberrin on 19/8/20.
//


#include "kgl_analysis_mutation_synthetic.h"
#include "kgl_analysis_mutation_syngen.h"
#include "kgl_analysis_mutation_inbreed_output.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


// Calculate the Synthetic Population Inbreeding Coefficient
bool kgl::SyntheticAnalysis::syntheticInbreeding(  std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   InbreedingParameters& parameters) {


  static const ContigOffset_t sampling_distance = 10000;
  // Filter out any variants that did not pass VCF filters (otherwise we get duplicate variants).
  unphased_ptr = unphased_ptr->filterVariants(PassFilter());

  parameters.lociiArguments().lociiSpacing(sampling_distance);

  syntheticInbreedingSample( unphased_ptr, parameters);


  return true;

}


bool kgl::SyntheticAnalysis::syntheticInbreedingSample( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                        InbreedingParameters& parameters) {


  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr, parameters.lociiArguments());

  ResultsMap results_map = processSynResults(contig_locus_map, parameters);

  InbreedingOutput::writeSynResults(results_map, parameters);

  return true;

}



kgl::ResultsMap kgl::SyntheticAnalysis::processSynResults( const ContigLocusMap& contig_locus_map,
                                                           InbreedingParameters& parameters) {

  ResultsMap results_map;
  // Retrieve the algorithm function object.
  auto algorithm_opt = InbreedingCalculation::namedAlgorithm(parameters.inbreedingAlgorthim());

  if (not algorithm_opt) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Inbreeding algorithm not found: {}", parameters.inbreedingAlgorthim());
    return  results_map;

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
      for (auto &future : future_vector) {

        auto locus_results = future.get();
        ExecEnv::log().info("Processed Diploid genome: {} for inbreeding and relatedness", locus_results.genome);
        results_map[locus_results.genome] = locus_results;

      }

    } // for locus

  } // for contig.

  return results_map;

}

