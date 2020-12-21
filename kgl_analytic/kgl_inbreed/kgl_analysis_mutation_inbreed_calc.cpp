//
// Created by kellerberrin on 21/8/20.
//

#include "kgl_analysis_mutation_inbreed_locus.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kel_distribution.h"

#include <fstream>
#include <algorithm>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


bool kgl::RetryCalcResult::checkRetry(double retry) {

  ++retry_count_;
  if (retry_count_ > max_retry_) {

    return true;

  }

  current_retries_.push_back(retry);

  if (current_retries_.size() > min_retry_) {

    current_retries_.pop_front();
    return checkTolerance();

  } else if (current_retries_.size() == min_retry_) {

    return checkTolerance();

  }

  return false;

}


bool kgl::RetryCalcResult::checkTolerance() const {

  auto current_entry = current_retries_.begin();
  while (current_entry != current_retries_.end()) {

    auto next_entry = ++current_entry;
    if (next_entry == current_retries_.end()) {

      break;

    }

    if (std::fabs(*current_entry - *next_entry) > tolerance_) {

      return false;

    }

    current_entry = next_entry;

  }

  return true;

}


std::optional<kgl::InbreedingAlgorithm> kgl::InbreedingCalculation::namedAlgorithm(const std::string& algorithm_name) {

  auto inbreeding_algo = inbreeding_algo_map_.find(algorithm_name);

  if (inbreeding_algo != inbreeding_algo_map_.end()) {

    auto [algo_name, algo_function] = *inbreeding_algo;
    return algo_function;

  }

  return std::nullopt;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These functions calculate the inbreeding coefficient for an individual by looking at multiple locii.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double kgl::InbreedingCalculation::logLikelihood(std::vector<double>& x, std::vector<AlleleFreqInfo>& data) {

  double log_prob_sum{0.0};
  double f = x[0];
  static const double small_prob{1e-10};

  for (auto const& allele_freq : data) {

    switch(allele_freq.alleleType()) {

      case AlleleClassType::MAJOR_HOMOZYGOUS:
      case AlleleClassType::MINOR_HOMOZYGOUS: {
        double freq_sqd = allele_freq.firstAllele().frequency() * allele_freq.firstAllele().frequency();
        double prob = (f * allele_freq.firstAllele().frequency()) + ((1.0 - f) * freq_sqd);
        prob = std::clamp<double>(prob, small_prob, 1.0);
        log_prob_sum += std::log(prob);
      }
        break;

      case AlleleClassType::MINOR_HETEROZYGOUS:
      case AlleleClassType::MAJOR_HETEROZYGOUS: {

        double prob = 2 * (1.0 - f) * allele_freq.firstAllele().frequency() * allele_freq.secondAllele().frequency();
        prob = std::clamp<double>(prob, small_prob, 1.0);
        log_prob_sum += std::log(prob);

      }
        break;

    }

  }

  return log_prob_sum;

}

kel::Optimize kgl::InbreedingCalculation::createLogLikelihoodOptimizer() {

  const size_t parameter_dimension = 1;
  Optimize optimizer(OptimizationAlgorithm::LN_NELDERMEAD, parameter_dimension, OptimizationType::MAXIMIZE);
  std::vector<double> lower_bound{-1.0};
  std::vector<double> upper_bound{1.0};
  optimizer.boundingHypercube(upper_bound, lower_bound);
  optimizer.stoppingCriteria(OptimizeStoppingType::ABSOLUTE_PARAMETER_THRESHOLD, {1e-06});
  optimizer.stoppingCriteria(OptimizeStoppingType::MAXIMUM_EVALUATIONS, { 500 });

  return optimizer;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The loglikelihood inbreeding algorithm.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::LocusResults
kgl::InbreedingCalculation::processLogLikelihood(const GenomeId_t& genome_id,
                                                 const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                 const std::string& super_population_field,
                                                 const std::shared_ptr<const ContigVariant>& locus_list,
                                                 const InbreedingParameters& parameters) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;
  // The real  distribution [upper, lower] used to as a start point for the loglike algorithm.
  UniformRealDistribution initialize_distribution(INIT_UPPER_, INIT_LOWER_);
  // Get the locus frequencies.
  auto [frequency_vector, locus_results] = generateFrequencies(genome_id,
                                                               contig_ptr,
                                                               super_population_field,
                                                               locus_list,
                                                               parameters.lociiArguments().frequencySource());

  double updated_coefficient{0.0};
  Optimize likelihood_optimizer = createLogLikelihoodOptimizer();
  RetryCalcResult retry_results(FINAL_ACCURACY_, MIN_RETRIES_, MAX_RETRIES_);

  // The calculation loop using the non-linear optimizer.
  do {

    // Random start on the unit interval.
    std::vector<double> coefficient;
    double initial_f = initialize_distribution.random(entropy_mt.generator());
    coefficient.push_back(initial_f);

    auto [result_code, value, iterations] = likelihood_optimizer.optimize<std::vector<AlleleFreqInfo>>( coefficient,
                                                                                                        frequency_vector,
                                                                                                        &logLikelihood);

    if (not Optimize::returnSuccess(result_code)) {

      ExecEnv::log().error("InbreedingAnalysis::processLogLikelihood; Genome: {}, loglikelihood: {}, Max at kgl_inbreed: {}, initial kgl_inbreed: {}, calc kgl_inbreed: {}, optimizer result: {}, iterations: {}",
                          locus_results.genome, value, coefficient.front(), initial_f, coefficient.front(), Optimize::returnDescription(result_code), iterations);

    }

    updated_coefficient = coefficient.front();

  } while (not retry_results.checkRetry(updated_coefficient));

  if (retry_results.retries() >= MAX_RETRIES_) {

    ExecEnv::log().warn( "InbreedingCalculation::processLogLikelihood, Genome: {}, retries: {}, final value: {} Loglikelihood Inbreeding algorithm did not converge",
                         locus_results.genome, retry_results.retries(), updated_coefficient);
    updated_coefficient = 0.0;

  }

  locus_results.inbred_allele_sum = updated_coefficient;

  ExecEnv::log().info("LogLikelihood: Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, IBD Inbreeding: {}, retries: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.minor_homo_count, locus_results.total_allele_count,
                      locus_results.inbred_allele_sum, retry_results.retries());

  return locus_results;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the Hall maximum expectation algorithm
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LocusResults
kgl::InbreedingCalculation::processHallME( const GenomeId_t& genome_id,
                                           const std::shared_ptr<const DiploidContig>& contig_ptr,
                                           const std::string& super_population_field,
                                           const std::shared_ptr<const ContigVariant>& locus_list,
                                           const InbreedingParameters& parameters) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;
  // The real unit distribution [upper, 0] used to as a start point for the EM algorithm.
  UniformRealDistribution initialize_distribution(INIT_UPPER_, 0);
  // Get the locus frequencies.
  auto [frequency_vector, locus_results] = generateFrequencies(genome_id,
                                                               contig_ptr,
                                                               super_population_field,
                                                               locus_list,
                                                               parameters.lociiArguments().frequencySource());

  double updated_coefficient{0.0};
  double inbreed_coefficient;
  RetryCalcResult retry_results(FINAL_ACCURACY_, MIN_RETRIES_, MAX_RETRIES_);

  do {

    // Random start on the unit interval.
    updated_coefficient = initialize_distribution.random(entropy_mt.generator());
    RetryCalcResult converge_retry(FINAL_ACCURACY_, MINIMUM_ITERATIONS_, MAXIMUM_ITERATIONS_);

    // Perform the EM algorithm
    do {

      inbreed_coefficient = updated_coefficient;
      double expectation_sum = 0.0;

      for (auto const& allele_freq : frequency_vector) {

        switch(allele_freq.alleleType()) {

          case AlleleClassType::MAJOR_HOMOZYGOUS:
          case AlleleClassType::MINOR_HOMOZYGOUS: {

            double denominator = (inbreed_coefficient + ((1.0-inbreed_coefficient) * allele_freq.firstAllele().frequency()));
            if (denominator != 0) {

              expectation_sum += inbreed_coefficient / denominator;

            }

          }
            break;

          case AlleleClassType::MAJOR_HETEROZYGOUS:
          case AlleleClassType::MINOR_HETEROZYGOUS:
            break;

        }

      }

      updated_coefficient = expectation_sum / static_cast<double>(frequency_vector.size());

    } while(not converge_retry.checkRetry(updated_coefficient));

  } while (not retry_results.checkRetry(updated_coefficient)); // while retries.

  if (retry_results.retries() >= MAX_RETRIES_) {

    ExecEnv::log().warn( "InbreedingCalculation::processHallME, Genome: {}, retries: {} Hall EM Inbreeding algorithm did not converge",
                         locus_results.genome, retry_results.retries());
    updated_coefficient = 0.0;

  }

  locus_results.inbred_allele_sum = updated_coefficient;

  ExecEnv::log().info("HallMe: Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, IBD Inbreeding: {}, retries: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.minor_homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum, retry_results.retries());

  return locus_results;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is multi locus inbreeding calculation based on the inbreeding definition (F = 1 - (Observed_Hetero/Expected_Hetero).
// Using the identity Expected_Hetero = TotalAlleles - Expected_Homo, this can be calculated using Expected_Homo or Expected_Hetero.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LocusResults
kgl::InbreedingCalculation::processSimple(const GenomeId_t& genome_id,
                                          const std::shared_ptr<const DiploidContig>& contig_ptr,
                                          const std::string& super_population_field,
                                          const std::shared_ptr<const ContigVariant>& locus_list,
                                          const InbreedingParameters& parameters) {

  // Get the locus frequencies.
  auto [frequency_vector, locus_results] = generateFrequencies( genome_id,
                                                                contig_ptr,
                                                                super_population_field,
                                                                locus_list,
                                                                parameters.lociiArguments().frequencySource());
  const bool calc_hetero{false}; // Which calculation to use.
  double heterozygous_inbreeding{0.0};
  double homozygous_inbreeding{0.0};

  if (locus_results.total_allele_count > 0) {

    auto observed_heterozygous = static_cast<double>(locus_results.major_hetero_count + locus_results.minor_hetero_count);
    auto observed_homozygous = static_cast<double>(locus_results.minor_homo_count + locus_results.major_homo_count);
    auto expected_heterozygous = locus_results.minor_hetero_freq + locus_results.major_hetero_freq;
    auto expected_homozygous = locus_results.minor_homo_freq + locus_results.major_homo_freq;

      // Heterozygous inbreeding
    heterozygous_inbreeding = 1.0 - (observed_heterozygous / expected_heterozygous);

      // Homozygous inbreeding
    homozygous_inbreeding = (observed_homozygous - expected_homozygous) / (static_cast<double>(locus_results.total_allele_count) - expected_homozygous);

    ExecEnv::log().info("Simple: Genome: {}, Exp Het: {}, Exp Hom: {}, Exp Het + Hom: {}, Obs Het: {}, Obs Hom: {}, Obs Het + Hom: {}, Het F: {}, Hom F: {}",
                        locus_results.genome, expected_heterozygous, expected_homozygous, (expected_heterozygous + expected_homozygous),
                        observed_heterozygous, observed_homozygous, locus_results.total_allele_count,
                        heterozygous_inbreeding, homozygous_inbreeding);

  }

  if (calc_hetero) {

    locus_results.inbred_allele_sum = heterozygous_inbreeding;

  } else {

    locus_results.inbred_allele_sum = homozygous_inbreeding;

  }

  return locus_results;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the published Ritland multi-locus calculation.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LocusResults
kgl::InbreedingCalculation::processRitlandLocus(const GenomeId_t &genome_id,
                                           const std::shared_ptr<const DiploidContig>& contig_ptr,
                                           const std::string& super_population_field,
                                           const std::shared_ptr<const ContigVariant>& locus_list,
                                           const InbreedingParameters& parameters) {

  // Ignore rare homozygous combinations, as these 'blow up' the ratio below.
  constexpr const static double minimum_frequency{0.001};

  size_t sum_allele{0};
  double locus_allele_sum{0.0};
  auto [frequency_vector, locus_results] = generateFrequencies( genome_id,
                                                                contig_ptr,
                                                                super_population_field,
                                                                locus_list,
                                                                parameters.lociiArguments().frequencySource());

  for (auto const& allele_freq : frequency_vector) {

    switch(allele_freq.alleleType()) {

      case AlleleClassType::MAJOR_HOMOZYGOUS:
      case AlleleClassType::MINOR_HOMOZYGOUS: {

        if (allele_freq.firstAllele().frequency() > minimum_frequency) {
          // This ratio becomes unstable for low frequency homozygous alleles.
          double ratio = (1.0 / allele_freq.firstAllele().frequency());
          locus_allele_sum += ratio;
          locus_allele_sum -= 1.0;
          ++sum_allele;

        }

      }
        break;

      case AlleleClassType::MAJOR_HETEROZYGOUS:
      case AlleleClassType::MINOR_HETEROZYGOUS: {

        locus_allele_sum -= 1.0;
        ++sum_allele;

      }
        break;

    }

  }


  locus_results.inbred_allele_sum = (sum_allele > 0 ? locus_allele_sum / static_cast<double>(sum_allele) : 0.0);

  ExecEnv::log().info("RitlandLocus: Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.minor_homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}


