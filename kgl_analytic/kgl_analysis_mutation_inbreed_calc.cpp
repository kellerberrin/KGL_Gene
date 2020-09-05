//
// Created by kellerberrin on 21/8/20.
//

#include "kel_thread_pool.h"
#include "kgl_analysis_mutation_inbreed_aux.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kel_distribution.h"

#include <fstream>

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These functions calculate the inbreeding coefficient for an individual by looking at multiple locii.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  locus_results.inbred_allele_sum = 0.0;
  locus_results.homo_count = 0;
  locus_results.hetero_count = 0;
  locus_results.total_allele_count = 0;
  std::vector<std::pair<double, double>> frequency_vector;

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
            ++locus_results.hetero_count;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset[0]->homozygous(*diploid_offset[1])) {


              auto[result, AF_value] = InbreedSampling::processFloatField(*locus_variant, super_population_field);
              if (result and AF_value > 0.0 and AF_value <= 1.0) {

                ++locus_results.homo_count;
                frequency_vector.emplace_back(AF_value, (AF_value * AF_value));

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

  double inbreed_coefficient;
  // Random start on the unit interval.
  double updated_coefficient = unit_distribution.random(entropy_mt.generator());
  size_t iteration_count = 0;

  ExecEnv::log().info( "Genome: {}, Frequency vector size: {}", locus_results.genome, frequency_vector.size());

  // Perform the EM algorithm
  do {

    inbreed_coefficient = updated_coefficient;
    double expectation_sum = 0.0;

    for (auto const& [frequency, sqr_frequency] : frequency_vector) {

      double fp = inbreed_coefficient * frequency;
      double denominator = (fp + ((1.0-inbreed_coefficient) * sqr_frequency));
      if (denominator != 0) {

        expectation_sum += fp / denominator;

      }

    }


    if (iteration_count > MAXIMUM_ITERATIONS_) {

      ExecEnv::log().warn( "InbreedingCalculation::processHallME, Genome: {}, Iterations: {} EM Inbreeding algorithm did not converge",
                           locus_results.genome, iteration_count);
      updated_coefficient = 0.0;
      break;

    }

    ++iteration_count;
    updated_coefficient = expectation_sum / static_cast<double>(locus_results.total_allele_count);

  } while(iteration_count < MINIMUM_ITERATIONS_ or std::fabs(updated_coefficient - inbreed_coefficient) > FINAL_ACCURACY_);

  locus_results.inbred_allele_sum = updated_coefficient;

  ExecEnv::log().info("Genome: {}, Super: {}, Het: {}, Hom: {}, Allele Count: {}, Inbreeding: {}, iterations: {}",
                      locus_results.genome, super_population_field, locus_results.hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum, iteration_count);

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
  locus_results.hetero_count = 0;
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
            ++locus_results.hetero_count;

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
                      locus_results.genome, super_population_field, locus_results.hetero_count,
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
        ++locus_results.hetero_count;

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
                      locus_results.genome, super_population_field, locus_results.hetero_count,
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
        ++locus_results.hetero_count;

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
                      locus_results.genome, super_population_field, locus_results.hetero_count,
                      locus_results.homo_count, locus_results.total_allele_count, locus_results.inbred_allele_sum);

  return locus_results;

}

