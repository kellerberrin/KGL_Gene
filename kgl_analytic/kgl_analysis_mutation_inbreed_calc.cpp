//
// Created by kellerberrin on 21/8/20.
//

#include "kgl_analysis_mutation_inbreed_aux.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kel_distribution.h"

#include <fstream>
#include <algorithm>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These functions calculate the inbreeding coefficient for an individual by looking at multiple locii.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double kgl::InbreedingCalculation::logLikelihood(std::vector<double>& x, std::vector<AlleleFreqInfo>& data) {

  double log_prob_sum = 0.0;
  double f = x[0];
  const double small_prob = 1e-10;
  for (auto const& [allele_type, allele1_prob, allele2_prob] : data) {

    switch(allele_type) {

      case MinorAlleleType::HOMOZYGOUS: {
        double prob = (f * allele1_prob) + ((1.0 - f) * allele1_prob * allele1_prob);
        prob = std::clamp<double>(prob, small_prob, 1.0);
        log_prob_sum += std::log(prob);
      }
        break;

      case MinorAlleleType::HETEROZYGOUS: {

        double prob = 2 * (1.0 - f) * allele1_prob * allele2_prob;
        prob = std::clamp<double>(prob, small_prob, 1.0);
        log_prob_sum += std::log(prob);

      }
        break;

      case MinorAlleleType::MINOR_HETEROZYGOUS: {
        double prob = 2 * (1.0 - f) * allele1_prob * allele2_prob;
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
// Sample locii data.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<kgl::AlleleFreqInfo>
kgl::InbreedingCalculation::generateGnomadFreq(const GenomeId_t& genome_id,
                                               const std::shared_ptr<const DiploidContig>& contig_ptr,
                                               const std::string& super_population_field,
                                               const std::shared_ptr<const ContigVariant>& locus_list) {

  std::vector<AlleleFreqInfo> frequency_vector;

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // Join on the diploid contig.
    auto locus_variant_array = offset_ptr->getVariantArray();
    // Check the diploid genome for any minor alleles at this location.
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);
    // If minor alleles are at the location.
    if (diploid_variant_opt) {

      auto const &diploid_offset = diploid_variant_opt.value();
      // Loop through the variants in the locus..
      for (auto const &locus_variant : locus_variant_array) {

        // The diploid offset can have 1 or 2 minor alleles
        // We need only check the first against the locus which is assumed to list all possible minor alleles.
        if (diploid_offset.front()->analogous(*locus_variant)) {
          // Found the matching locus allele.

          auto[result, AF_value] = InbreedSampling::processFloatField(*locus_variant, super_population_field);
          if (not result) {

            // Problem obtaining allele frequency.
            ExecEnv::log().warn("InbreedingAnalysis::processLogLikelihood; Genome: {}, Unable to obtain Locus allele frequency: {} for SNP: {}",
                                genome_id, super_population_field, locus_variant->output(',', VariantOutputIndex::START_0_BASED, false));

          } else if (AF_value <= 0.0 or AF_value >= 1.0) {

            continue;

          } else { // We have a valid minor allele frequency

            if (diploid_offset.size() == 1) {
              // The sample is alt allele heterozygous
              // Calculate the Major allele frequency as the complement of the sum of minor alleles.
              double minor_allele_frequency_sum = 0;
              for (auto const& minor_variant : locus_variant_array) {

                auto[minor_result, minor_frequency] = InbreedSampling::processFloatField(*minor_variant, super_population_field);
                if (not minor_result) {

                  // Problem obtaining allele frequency.
                  ExecEnv::log().warn("InbreedingAnalysis::processLogLikelihood; Genome: {}, Unable to obtain Locus allele frequency: {} for SNP: {}",
                                      genome_id, super_population_field, minor_variant->output(',', VariantOutputIndex::START_0_BASED, false));

                } else {


                  minor_allele_frequency_sum += minor_frequency;

                }

              }

              double major_allele_frequency = std::clamp<double>((1.0 - minor_allele_frequency_sum), 0.0, 1.0);
              // Add to the loglikelihood vector.
              frequency_vector.push_back({MinorAlleleType::HETEROZYGOUS, AF_value, major_allele_frequency});

            } else if (diploid_offset.size() == 2) {

              if (diploid_offset[0]->homozygous(*diploid_offset[1])) {

                // Add to the loglikelihood vector.
                frequency_vector.push_back({MinorAlleleType::HOMOZYGOUS, AF_value, AF_value});


              } else {
                // The sample has different alt alleles.
                // Find the minor allele and obtain it's frequency.
                bool found_second_minor = false;
                for (auto& second_variant : locus_variant_array) {

                  if (diploid_offset[1]->analogous(*second_variant)) {

                    // Found the second variant.
                    auto[minor_result, minor_frequency] = InbreedSampling::processFloatField(*second_variant, super_population_field);
                    if (not minor_result) {

                      // Problem obtaining allele frequency.
                      ExecEnv::log().warn("InbreedingAnalysis::processLogLikelihood; Genome: {}, Unable to obtain Locus allele frequency: {} for SNP: {}",
                                          genome_id, super_population_field, second_variant->output(',', VariantOutputIndex::START_0_BASED, false));

                    } else {

                      minor_frequency = std::clamp(minor_frequency, 0.0, 1.0);
                      frequency_vector.push_back({MinorAlleleType::MINOR_HETEROZYGOUS, AF_value, minor_frequency});
                      found_second_minor = true;

                    }

                    break;

                  }

                }

                if (not found_second_minor) {
                  ExecEnv::log().warn("InbreedingCalculation::processLogLikelihood; Genome: {}, Not Found Second Minor SNP: {}",
                                      genome_id, diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false));

                }

              } // Two minor alleles

            }  // if 2 alleles

            break; // No need to search further.

          } // Valid minor allele frequency.

        } // Found locus variant

      } // For all locus variants

    } // found matching allele.

  } // for all locii

  return frequency_vector;

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
                                                 const std::shared_ptr<const ContigVariant>& locus_list) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;
  // The real unit distribution [1,0] used to as a start point for the loglik algorithm.
  UniformUnitDistribution unit_distribution;

  LocusResults locus_results;
  locus_results.genome = genome_id;

  auto frequency_vector = generateGnomadFreq(genome_id, contig_ptr, super_population_field, locus_list);

  double sum = 0.0;
  double sq_sum = 0.0;
  for (auto [homozygous, allele1_prob, allele2_prob] : frequency_vector) {

    ++locus_results.total_allele_count;

    switch(homozygous) {

      case MinorAlleleType::HOMOZYGOUS:
        ++locus_results.homo_count;
        break;

      case MinorAlleleType::HETEROZYGOUS:
        ++locus_results.major_hetero_count;
        break;

      case MinorAlleleType::MINOR_HETEROZYGOUS:
        ++locus_results.minor_hetero_count;
        break;

    }

    sum += allele1_prob;
    sq_sum += (allele1_prob * allele1_prob);

  }
  double p = sum / static_cast<double>(frequency_vector.size());
  double binomial_mean = static_cast<double>(frequency_vector.size()) * p;
  double binomial_var = binomial_mean * (1.0 - p);
  double poisson_binomial_mean = sum;
  double poisson_binomial_var = sum - sq_sum;

  ExecEnv::log().info( "Genome: {}, Frequency vector size: {}, mean prob: {}, bin mean: {} var: {}, poisson bin mean: {}, var: {}",
                       locus_results.genome, frequency_vector.size(), p, binomial_mean, binomial_var, poisson_binomial_mean, poisson_binomial_var);


  double updated_coefficient = 0.0;
  double inbreed_coefficient;
  size_t retries = 0;
  Optimize likelihood_optimizer = createLogLikelihoodOptimizer();

  do {


    // Random start on the unit interval.
    std::vector<double> coefficient;
    double initial_f = unit_distribution.random(entropy_mt.generator());
    coefficient.push_back(initial_f);

    inbreed_coefficient = updated_coefficient;
    ++retries;

    auto [result_code, value, iterations] = likelihood_optimizer.optimize<std::vector<AlleleFreqInfo>>( coefficient,
                                                                                                        frequency_vector,
                                                                                                        &logLikelihood);

    if (not Optimize::returnSuccess(result_code)) {

      ExecEnv::log().error("InbreedingAnalysis::processLogLikelihood; Genome: {}, loglikelihood: {}, Max at inbreed: {}, initial inbreed: {}, previous inbreed: {}, optimizer result: {}, iterations: {}",
                          locus_results.genome, value, coefficient.front(), initial_f, inbreed_coefficient, Optimize::returnDescription(result_code), iterations);

    }

    updated_coefficient = coefficient.front();

  } while (retries < MIN_RETRIES_ or (retries <= MAX_RETRIES_ and std::fabs(updated_coefficient - inbreed_coefficient) >= 1e-04));

  if (retries >= MAX_RETRIES_) {

    ExecEnv::log().warn( "InbreedingCalculation::processLogLikelihood, Genome: {}, retries: {}, final value: {} Loglikelihood Inbreeding algorithm did not converge",
                         locus_results.genome, retries, updated_coefficient);
    updated_coefficient = 0.0;

  }

  locus_results.inbred_allele_sum = updated_coefficient;

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, IBD Inbreeding: {}, retries: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum, retries);


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
                                           const std::shared_ptr<const ContigVariant>& locus_list) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;
  // The real unit distribution [1,0] used to as a start point for the EM algorithm.
  UniformUnitDistribution unit_distribution;

  LocusResults locus_results;
  locus_results.genome = genome_id;

  auto frequency_vector = generateGnomadFreq(genome_id, contig_ptr, super_population_field, locus_list);

  double sum = 0.0;
  double sq_sum = 0.0;
  for (auto [homozygous, allele1_prob, allele2_prob] : frequency_vector) {

    ++locus_results.total_allele_count;

    switch(homozygous) {

      case MinorAlleleType::HOMOZYGOUS:
        ++locus_results.homo_count;
        break;

      case MinorAlleleType::HETEROZYGOUS:
        ++locus_results.major_hetero_count;
        break;

      case MinorAlleleType::MINOR_HETEROZYGOUS:
        ++locus_results.minor_hetero_count;
        break;

    }

    sum += allele1_prob;
    sq_sum += (allele1_prob * allele1_prob);

  }


  double p = sum/static_cast<double>(frequency_vector.size());
  double binomial_mean = static_cast<double>(frequency_vector.size()) * p;
  double binomial_var = binomial_mean * (1.0 - p);
  double poisson_binomial_mean = sum;
  double poisson_binomial_var = sum - sq_sum;


  ExecEnv::log().info( "Genome: {}, Frequency vector size: {}, mean prob: {}, bin mean: {} var: {}, poisson bin mean: {}, var: {}",
                       locus_results.genome, frequency_vector.size(), p, binomial_mean, binomial_var, poisson_binomial_mean, poisson_binomial_var);

  double updated_coefficient = 0.0;
  double inbreed_coefficient;
  double previous_coefficient;
  size_t retries = 0;

  do {

    // Random start on the unit interval.
    previous_coefficient = updated_coefficient;
    updated_coefficient = unit_distribution.random(entropy_mt.generator());

    size_t iteration_count = 0;
    // Perform the EM algorithm
    do {

      inbreed_coefficient = updated_coefficient;
      double expectation_sum = 0.0;

      for (auto [homozygous, allele1_prob, allele2_prob] : frequency_vector) {

        switch(homozygous) {

          case MinorAlleleType::HOMOZYGOUS: {

            double denominator = (inbreed_coefficient + ((1.0-inbreed_coefficient) * allele1_prob));
            if (denominator != 0) {

              expectation_sum += inbreed_coefficient / denominator;

            }

          }
            break;

          case MinorAlleleType::HETEROZYGOUS:
          case MinorAlleleType::MINOR_HETEROZYGOUS:
            break;

        }

      }

      ++iteration_count;
      updated_coefficient = expectation_sum / static_cast<double>(frequency_vector.size());

    } while(iteration_count < MINIMUM_ITERATIONS_
            or (iteration_count < MAXIMUM_ITERATIONS_
            and std::fabs(updated_coefficient - inbreed_coefficient) > FINAL_ACCURACY_));

    ++retries;

  } while (retries < MIN_RETRIES_
           or (retries < MAX_RETRIES_
           and std::fabs(updated_coefficient - previous_coefficient) > FINAL_ACCURACY_)); // while retries.

  if (retries >= MAX_RETRIES_) {

    ExecEnv::log().warn( "InbreedingCalculation::processHallME, Genome: {}, retries: {} Hall EM Inbreeding algorithm did not converge",
                         locus_results.genome, retries);
    updated_coefficient = 0.0;

  }

  locus_results.inbred_allele_sum = updated_coefficient;

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, IBD Inbreeding: {}, retries: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum, retries);

  return locus_results;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the experimental Ritland multi-locus calculation developed in my document.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LocusResults
kgl::InbreedingCalculation::processRitlandMME( const GenomeId_t& genome_id,
                                               const std::shared_ptr<const DiploidContig>& contig_ptr,
                                               const std::string& super_population_field,
                                               const std::shared_ptr<const ContigVariant>& locus_list) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  LocusResults locus_results;
  locus_results.genome = genome_id;
  locus_results.inbred_allele_sum = 0.0;
  locus_results.homo_count = 0;
  locus_results.major_hetero_count = 0;
  locus_results.total_allele_count = 0;

  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // Join on the diploid contig.
    auto locus_variant_array = offset_ptr->getVariantArray();
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {

      auto const &diploid_offset = diploid_variant_opt.value();

      for (auto const &locus_variant : locus_variant_array) {

        if (diploid_offset[0]->analogous(*locus_variant)) {
          // Found the matching locus allele.
          // Get the allele super population frequency
          ++locus_results.total_allele_count;

          if (diploid_offset.size() == 1) {
            // The sample is alt allele heterozygous
            ++locus_results.major_hetero_count;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset[0]->homozygous(*diploid_offset[1])) {


              auto[result, AF_value] = InbreedSampling::processFloatField(*locus_variant, super_population_field);
              if (result and AF_value >= 0.0 and AF_value <= 1.0) {

                locus_results.inbred_allele_sum += (AF_value * AF_value);
                ++locus_results.homo_count;

              }

            } else {
              // The sample has different alt alleles. Possible but unlikely.
              ExecEnv::log().info("InbreedingAnalysis::multiLocus1; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                                  genome_id,
                                  diploid_offset[0]->output(',', VariantOutputIndex::START_0_BASED, false),
                                  diploid_offset[1]->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          }  // if 2 alleles

          break; // No need to search further.

        } // Found locus variant

      } // for all locus variants

    } // found matching allele.

  } // for all locii

  if (locus_results.total_allele_count > 0) {

    double homo_average = static_cast<double>(locus_results.homo_count) / static_cast<double>(locus_results.total_allele_count);
    double expected_average = locus_results.inbred_allele_sum /  static_cast<double>(locus_results.total_allele_count);
    locus_results.inbred_allele_sum = (homo_average - expected_average) / (1 - expected_average);

  } else {

    locus_results.inbred_allele_sum = 0.0;

  }

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

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
                                           const std::shared_ptr<const ContigVariant>& locus_list) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  LocusResults locus_results;
  locus_results.inbred_allele_sum = 0.0;
  locus_results.genome = genome_id;
  size_t sum_alternate_allele = 0;
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {
    // Join on the diploid contig.

    auto locus_variant_array = offset_ptr->getVariantArray();
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {

      sum_alternate_allele += locus_variant_array.size();
      ++locus_results.total_allele_count;
      // Determine if the sample alternate allele is Hom/Het or Mixed.
      auto const& diploid_offset = diploid_variant_opt.value();
      if (diploid_offset.size() == 1) {
        // The sample is alt allele heterozygous
        ++locus_results.major_hetero_count;

      } else if (diploid_offset.size() == 2) {

         if (diploid_offset[0]->homozygous(*diploid_offset[1])) {
          // The sample is alt allele homozygous
          // Find the matching locus allele
           ++locus_results.homo_count;

          for (auto const& locus_variant : locus_variant_array) {

            if (diploid_offset[0]->analogous(*locus_variant)) {
              // Found the matching locus allele.
              // Get the allele super population frequency
              auto [result, AF_value] = InbreedSampling::processFloatField(*locus_variant, super_population_field);
              if (result and AF_value > 0.0 and AF_value < 1.0) {

                locus_results.inbred_allele_sum += (1.0 / AF_value);
                locus_results.inbred_allele_sum -= 1.0;

              } // valid AF

              break; // No need to search further.

            } // Found locus variant

          } // for all locus variants

        } else {
          // The sample has different alt alleles.
          ExecEnv::log().info("InbreedingAnalysis::processRitlandLocus; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                              genome_id,
                              diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false),
                              diploid_offset[1]->output(',',VariantOutputIndex::START_0_BASED, false));

        }

      } else {

        ExecEnv::log().error("InbreedingAnalysis::processRitlandLocus; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
                             genome_id, diploid_offset.size(), offset, contig_ptr->contigId());
        continue;
      }

    }

  } // for all locus variants.

  locus_results.inbred_allele_sum = (sum_alternate_allele > 0 ? locus_results.inbred_allele_sum / static_cast<double>(sum_alternate_allele) : 0.0);

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the published Ritland multi-locus calculation.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::LocusResults
kgl::InbreedingCalculation::processRitlandPopulation( const GenomeId_t &genome_id,
                                                   const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                   const std::string& super_population_field,
                                                   const std::shared_ptr<const ContigVariant>& locus_list) {

  // Only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  LocusResults locus_results;
  locus_results.inbred_allele_sum = 0.0;
  locus_results.genome = genome_id;
  size_t sum_alternate_allele = 0;
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {
    // Join on the diploid contig.

    auto locus_variant_array = offset_ptr->getVariantArray();
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);

    if (diploid_variant_opt) {

      sum_alternate_allele += locus_variant_array.size();
      ++locus_results.total_allele_count;
      // Determine if the sample alternate allele is Hom/Het or Mixed.
      auto const& diploid_offset = diploid_variant_opt.value();
      if (diploid_offset.size() == 1) {
        // The sample is alt allele heterozygous
        ++locus_results.major_hetero_count;

      } else if (diploid_offset.size() == 2) {

        if (diploid_offset[0]->homozygous(*diploid_offset[1])) {
          // The sample is alt allele homozygous
          // Find the matching locus allele
          ++locus_results.homo_count;

          for (auto const& locus_variant : locus_variant_array) {

            if (diploid_offset[0]->analogous(*locus_variant)) {
              // Found the matching locus allele.
              // Get the allele super population frequency
              auto [result, AF_value] = InbreedSampling::processFloatField(*diploid_offset[0],
                                                                           InbreedSampling::inverseSuperPopulationField(super_population_field));
              if (result and AF_value > 0.0 and AF_value < 1.0) {

                locus_results.inbred_allele_sum += (1.0 / AF_value);
                locus_results.inbred_allele_sum -= 1.0;

              } // valid AF

              break; // No need to search further.

            } // Found locus variant

          } // for all locus variants

        } else {
          // The sample has different alt alleles.
          ExecEnv::log().info("InbreedingAnalysis::processRitlandLocus; Diploid genome: {} has two different non-ref alleles\n{}\n{}",
                              genome_id,
                              diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false),
                              diploid_offset[1]->output(',',VariantOutputIndex::START_0_BASED, false));

        }

      } else {

        ExecEnv::log().error("InbreedingAnalysis::processRitlandLocus; Diploid genome: {} has: {} SNPs at offset: {} contig: {}",
                             genome_id, diploid_offset.size(), offset, contig_ptr->contigId());
        continue;
      }

    }

  } // for all locus variants.

  locus_results.inbred_allele_sum = (sum_alternate_allele > 0 ? locus_results.inbred_allele_sum / static_cast<double>(sum_alternate_allele) : 0.0);

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}",
                      locus_results.genome, super_population_field, locus_results.major_hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}

