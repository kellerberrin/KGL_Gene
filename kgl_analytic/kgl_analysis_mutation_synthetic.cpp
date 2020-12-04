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
                                                   const InbreedingParameters& parameters,
                                                   InbreedingOutputResults& results) {


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

  InbreedingParameters local_params = parameters;
  std::vector<ContigOffset_t> locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                                                VariantDatabaseRead::SUPER_POP_ALL_,
                                                                                local_params.lociiArguments());
  local_params.lociiArguments().upperOffset(locii_vector.back());

  while (local_params.lociiArguments().upperOffset() < parameters.lociiArguments().upperOffset() and locii_vector.size() >= 100) {

    ResultsMap results_map = processSynResults( unphased_ptr, parameters);

    results.insertResults({local_params, results_map});

    local_params.lociiArguments().lowerOffset(local_params.lociiArguments().upperOffset());
    locii_vector = RetrieveLociiVector::getLociiCount(contig_ptr,
                                                      VariantDatabaseRead::SUPER_POP_ALL_,
                                                      local_params.lociiArguments());
    local_params.lociiArguments().upperOffset(locii_vector.back());

  }

  return true;

}


kgl::ResultsMap kgl::SyntheticAnalysis::processSynResults( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                           const InbreedingParameters& parameters) {

  ResultsMap results_map;
  // Retrieve the algorithm function object.
  auto algorithm_opt = InbreedingCalculation::namedAlgorithm(parameters.inbreedingAlgorthim());

  if (not algorithm_opt) {

    ExecEnv::log().error("InbreedingAnalysis::populationInbreeding, Inbreeding algorithm not found: {}", parameters.inbreedingAlgorthim());
    return  results_map;

  }

  ContigLocusMap contig_locus_map = InbreedSampling::getPopulationLocusMap(unphased_ptr, parameters.lociiArguments());

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

