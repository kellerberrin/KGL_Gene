//
// Created by kellerberrin on 2/11/20.
//

#include "kgl_analysis_mutation_inbreed_aux.h"
#include "kgl_filter.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kel_distribution.h"

#include <fstream>
#include <algorithm>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;



double kgl::AlleleFreqInfo::majorAlleleFrequency() const {

  double sum_allele_freq{0.0};
  for (auto const& allele_freq : allele_frequencies_) {

    sum_allele_freq += allele_freq.frequency();

  }

  return std::clamp((1.0 - sum_allele_freq), 0.0, 1.0);

}

// The probability of heterozygous alleles. This includes MINOR_HETEROZYGOUS (rare, two heterozygous minor alleles).
double kgl::AlleleFreqInfo::probAllHeterozygous() const {

  double sum_hetero_freq{0.0};
  size_t minor_allele_count = allele_frequencies_.size();

  for (size_t idx1 = 0; idx1 < minor_allele_count; ++idx1) {

    for (size_t idx2 = (idx1 + 1); idx2 < minor_allele_count; ++idx2) {

        sum_hetero_freq += 2.0 * allele_frequencies_[idx1].frequency() * allele_frequencies_[idx2].frequency();

    }

  }

  sum_hetero_freq += probHeterozygous();

  return sum_hetero_freq;

}

// The probability of major heterozygous alleles. This is the case where only a single minor allele is recorded (common).
double kgl::AlleleFreqInfo::probHeterozygous() const {

  double sum_hetero_freq{0.0};
  double major_allele_freq = majorAlleleFrequency();
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_hetero_freq += 2.0 * minor_allele.frequency() * major_allele_freq;

  }

  return sum_hetero_freq;

}


// The probability of minor homozygous alleles
double kgl::AlleleFreqInfo::probMinorHomozygous() const {

  double sum_homozygous_freq{0.0};
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_homozygous_freq += minor_allele.frequency() * minor_allele.frequency();

  }

  return sum_homozygous_freq;

}


// The probability of minor and major homozygous alleles
double kgl::AlleleFreqInfo::probAllHomozygous() const {

  double sum_homozygous_freq{0.0};
  double sum_allele_freq{0.0};
  for (auto const& minor_allele  : allele_frequencies_) {

    double frequency = minor_allele.frequency();
    sum_allele_freq += frequency;
    sum_homozygous_freq += frequency * frequency;

  }

  // Probability of a homozygous major allele
  sum_homozygous_freq += (1.0 - sum_allele_freq) * (1.0 - sum_allele_freq);

  return sum_homozygous_freq;

}


double kgl::AlleleFreqInfo::sumFrequencies() const {

  double sum_freq{0.0};
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_freq += minor_allele.frequency();

  }

  return sum_freq;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate Allele Frequencies using Diploid.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::pair<std::vector<kgl::AlleleFreqInfo>, kgl::LocusResults>
kgl::InbreedingCalculation::generateDiploidFreq( const GenomeId_t& genome_id,
                                                 const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                 const std::string& super_population_field,
                                                 const std::shared_ptr<const ContigVariant>& locus_list) {

  std::vector<AlleleFreqInfo> frequency_vector;
  LocusResults locus_results;
  locus_results.genome = genome_id;

  // Convert to the super population field codes used in Gnomad.
  const std::string gnomad_super_pop = InbreedSampling::lookupSuperPopulationField(super_population_field);
  // The diploid super population field code.
  const std::string diploid_super_pop = super_population_field + "_AF";

  // Diploid contig, only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // Join on the diploid contig.
    auto locus_variant_array = offset_ptr->getVariantArray();
    // Check the diploid genome for any minor alleles at this location.
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);
    // If minor alleles are at the diploid location (rare if not).
    if (diploid_variant_opt) {

      auto const &diploid_offset = diploid_variant_opt.value();

      std::vector<AlleleFreqRecord> allele_vector;
      double minor_allele_frequencies{0.0};
      // Loop through the variants in the locus..
      for (auto const &locus_variant : locus_variant_array) {

        auto[result, AF_value] = InbreedSampling::processFloatField(*locus_variant, gnomad_super_pop);
        if (not result) {

          // Problem obtaining allele frequency.
          ExecEnv::log().warn("InbreedingCalculation::generateDiploidFreq; Genome: {}, Unable to obtain Gnomad allele frequency: {} for SNP: {}",
                              genome_id, gnomad_super_pop, locus_variant->output(',', VariantOutputIndex::START_0_BASED, false));

        } else  {

          AF_value = std::clamp(AF_value, 0.0, 1.0); // Just in case.
          // Store the Gnomad variant and super population frequency
          minor_allele_frequencies += AF_value;
          allele_vector.emplace_back(locus_variant, AF_value);

        }

      }

      // Loop through the reference (Gnomad) frequency vector
      for (auto const &allele_freq : allele_vector) {
        // The diploid offset can have 1 or 2 minor alleles
        // We need only check the first against the locus which is assumed to list all possible minor alleles.
        if (diploid_offset.front()->analogous(*allele_freq.allele())) {

          if (diploid_offset.size() == 1) {
            // The sample is minor allele heterozygous, read the allele frequency using the diploid variant.
            auto[result, AF_value] = InbreedSampling::processFloatField(*diploid_offset.front(), diploid_super_pop);
            if (not result) {

              // Problem obtaining allele frequency.
              ExecEnv::log().warn("InbreedingCalculation::generateDiploidFreq; Genome: {}, Unable to obtain Heterozygous Diploid allele frequency: {} for SNP: {}",
                                  genome_id, diploid_super_pop, diploid_offset.front()->output(',', VariantOutputIndex::START_0_BASED, false));

            } else  {

              AF_value = std::clamp(AF_value, 0.0, 1.0); // Just in case.
              AlleleFreqRecord minor_allele(diploid_offset.front(), AF_value);
              // Calculate the major allele frequency as the complement of the sum of minor alleles.
              // Adjusted for the diploid allele freq.
              minor_allele_frequencies -= allele_freq.frequency(); // Subtract the Gnomad frequency
              minor_allele_frequencies += AF_value; // Add in the diploid frequency.
              double major_allele_frequency = std::clamp<double>((1.0 - minor_allele_frequencies), 0.0, 1.0);
              AlleleFreqRecord second_allele(diploid_offset.front(), major_allele_frequency);
              // Add to the frequency vector.
              frequency_vector.emplace_back(MinorAlleleType::HETEROZYGOUS, minor_allele, second_allele, std::move(allele_vector));
              // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
            }
            break;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset[0]->homozygous(*diploid_offset[1])) {

              auto[result, AF_value] = InbreedSampling::processFloatField(*diploid_offset.front(), diploid_super_pop);
              if (not result) {

                // Problem obtaining allele frequency.
                ExecEnv::log().warn("InbreedingCalculation::generateDiploidFreq; Genome: {}, Unable to obtain Homozygous Diploid allele frequency: {} for SNP: {}",
                                    genome_id, diploid_super_pop, diploid_offset.front()->output(',', VariantOutputIndex::START_0_BASED, false));

              } else  {

                AF_value = std::clamp(AF_value, 0.0, 1.0); // Just in case.
                AlleleFreqRecord minor_allele(diploid_offset.front(), AF_value);
                AlleleFreqRecord second_allele(diploid_offset.back(), AF_value);
                frequency_vector.emplace_back(MinorAlleleType::HOMOZYGOUS, minor_allele, second_allele, std::move(allele_vector));
                // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
              }
              break;


            } else {
              // The sample has different alt alleles (rare).
              double minor_allele_freq{0.0};
              double second_allele_freq{0.0};

              auto[result, AF_value] = InbreedSampling::processFloatField(*diploid_offset.front(), diploid_super_pop);
              if (not result) {

                // Problem obtaining allele frequency.
                ExecEnv::log().warn("InbreedingCalculation::generateDiploidFreq; Genome: {}, Unable to obtain 1st Alt Diploid allele frequency: {} for SNP: {}",
                                    genome_id, diploid_super_pop, diploid_offset.front()->output(',', VariantOutputIndex::START_0_BASED, false));

                break;
              }
              minor_allele_freq = std::clamp(AF_value, 0.0, 1.0);

              auto[second_result, second_AF_value] = InbreedSampling::processFloatField(*diploid_offset.back(), diploid_super_pop);
              if (not second_result) {

                // Problem obtaining allele frequency.
                ExecEnv::log().warn("InbreedingCalculation::generateDiploidFreq; Genome: {}, Unable to 2nd Alt Diploid Locus allele frequency: {} for SNP: {}",
                                    genome_id, diploid_super_pop, diploid_offset.back()->output(',', VariantOutputIndex::START_0_BASED, false));

                break;
              }
              second_allele_freq = std::clamp(second_AF_value, 0.0, 1.0);

              AlleleFreqRecord minor_allele(diploid_offset.front(), minor_allele_freq);
              AlleleFreqRecord second_allele(diploid_offset.back(), second_allele_freq);
              frequency_vector.emplace_back(MinorAlleleType::MINOR_HETEROZYGOUS, minor_allele, second_allele, std::move(allele_vector));
              // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
              break;

            } // Two minor alleles

          }  // if 2 alleles

        } // Valid minor allele frequency.

      } // For all minor alleles in the reference at offset

    } // If minor diploid allele found at offset.

  } // For all offset locii.

  // Generate some frequency statistics.
  std::vector<double> freq_difference;
  for (auto const& allele_freq : frequency_vector) {

    ++locus_results.total_allele_count;

    switch(allele_freq.alleleType()) {

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

    bool found{false};
    for (auto const& ref_allele : allele_freq.alleleFrequencies()) {

      if (ref_allele.allele()->analogous(*allele_freq.minorAllele().allele())) {

        found = true;
        double difference = ref_allele.frequency() - allele_freq.minorAllele().frequency();
        freq_difference.push_back(difference);
        break;

      }

    } // End compare frequencies.

    // Should be a matching allele, complain if not.
    if (not found) {

      ExecEnv::log().error("InbreedingCalculation::generateDiploidFreq; Genome: {}, No matching reference SNP for Diploid SNP: {}",
                           genome_id, allele_freq.minorAllele().allele()->output(',', VariantOutputIndex::START_0_BASED, false));

    }

  }

  auto [mean, stddev] = Utility::stddev(freq_difference);
  ExecEnv::log().info( "Genome: {}, (Reference - Diploid) Allele frequency difference, Mean: {}, Stddev: {}", mean, stddev);

  return { frequency_vector, locus_results };

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate Allele Frequencies using Gnomad.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<kgl::AlleleFreqInfo>, kgl::LocusResults>
kgl::InbreedingCalculation::generateGnomadFreq(const GenomeId_t& genome_id,
                                               const std::shared_ptr<const DiploidContig>& contig_ptr,
                                               const std::string& super_population_field,
                                               const std::shared_ptr<const ContigVariant>& locus_list) {

  std::vector<AlleleFreqInfo> frequency_vector;
  LocusResults locus_results;
  locus_results.genome = genome_id;

  // Convert to the super population field codes used in Gnomad.
  const std::string gnomad_super_pop = InbreedSampling::lookupSuperPopulationField(super_population_field);

  // Diploid contig, only want SNP variants that have passed VCF filters.
  auto snp_contig_ptr = contig_ptr->filterVariants(AndFilter(SNPFilter(), PassFilter()));

  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // Join on the diploid contig.
    auto locus_variant_array = offset_ptr->getVariantArray();
    // Check the diploid genome for any minor alleles at this location.
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);
    // If minor alleles are at the diploid location (rare if not).
    if (diploid_variant_opt) {

      auto const &diploid_offset = diploid_variant_opt.value();

      std::vector<AlleleFreqRecord> allele_vector;
      double minor_allele_frequencies{0.0};
      // Loop through the variants in the locus..
      for (auto const &locus_variant : locus_variant_array) {

        auto[result, AF_value] = InbreedSampling::processFloatField(*locus_variant, gnomad_super_pop);
        if (not result) {

          // Problem obtaining allele frequency.
          ExecEnv::log().warn("InbreedingAnalysis::processLogLikelihood; Genome: {}, Unable to obtain Locus allele frequency: {} for SNP: {}",
                              genome_id, super_population_field, locus_variant->output(',', VariantOutputIndex::START_0_BASED, false));

        } else  {

          AF_value = std::clamp(AF_value, 0.0, 1.0); // Just in case.
          // Store the Gnomad variant and super population frequency
          minor_allele_frequencies += AF_value;
          allele_vector.emplace_back(locus_variant, AF_value);

        }

      }

      // Loop through the reference (Gnomad) frequency vector
      for (auto const &allele_freq : allele_vector) {
        // The diploid offset can have 1 or 2 minor alleles
        // We need only check the first against the locus which is assumed to list all possible minor alleles.
        if (diploid_offset.front()->analogous(*allele_freq.allele())) {

          if (diploid_offset.size() == 1) {
            // The sample is minor allele heterozygous
            AlleleFreqRecord minor_allele(diploid_offset.front(), allele_freq.frequency());
            // Calculate the major allele frequency as the complement of the sum of minor alleles.
            double major_allele_frequency = std::clamp<double>((1.0 - minor_allele_frequencies), 0.0, 1.0);
            AlleleFreqRecord second_allele(diploid_offset.front(), major_allele_frequency);
            // Add to the frequency vector.
            frequency_vector.emplace_back(MinorAlleleType::HETEROZYGOUS, minor_allele, second_allele, std::move(allele_vector));
            // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
            break;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset[0]->homozygous(*diploid_offset[1])) {

              AlleleFreqRecord minor_allele(diploid_offset.front(), allele_freq.frequency());
              AlleleFreqRecord second_allele(diploid_offset.back(), allele_freq.frequency());
              frequency_vector.emplace_back(MinorAlleleType::HOMOZYGOUS, minor_allele, second_allele, std::move(allele_vector));
              // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
              break;


            } else {
              // The sample has different alt alleles (rare).
              // Find the second minor allele and obtain it's frequency.
              bool found_second_minor = false;
              double second_frequency{0.0};
              for (auto const& second_allele_freq : allele_vector) {

                if (diploid_offset.back()->analogous(*second_allele_freq.allele())) {

                  second_frequency = second_allele_freq.frequency();
                  found_second_minor= true;
                  break;

                }

              } // For find second allele frequency.

              if (not found_second_minor) {
                ExecEnv::log().warn("InbreedingCalculation::processLogLikelihood; Genome: {}, Not Found Second Minor SNP: {}",
                                    genome_id, diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false));

              } else {

                AlleleFreqRecord minor_allele(diploid_offset.front(), allele_freq.frequency());
                AlleleFreqRecord second_allele(diploid_offset.back(), second_frequency);
                frequency_vector.emplace_back(MinorAlleleType::MINOR_HETEROZYGOUS, minor_allele, second_allele, std::move(allele_vector));
                // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
                break;

              }

            } // Two minor alleles

          }  // if 2 alleles

        } // Valid minor allele frequency.

      } // For all minor alleles in the reference at offset

    } // If minor diploid allele found at offset.

  } // For all offset locii.

  // Generate some frequency statistics.
  double sum{0.0};
  double sq_sum{0.0};

  locus_results.total_allele_count = frequency_vector.size();
  for (auto const& allele_freq : frequency_vector) {

    switch(allele_freq.alleleType()) {

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

    double freq = allele_freq.minorAllele().frequency();
    sum += freq;
    sq_sum += (freq * freq);

  }

  // Compare the first two moments assuming Binomial and Poisson Binomial distributions.
  double p = sum / static_cast<double>(frequency_vector.size());
  double binomial_mean = static_cast<double>(frequency_vector.size()) * p;
  double binomial_var = binomial_mean * (1.0 - p);
  double poisson_binomial_mean = sum;
  double poisson_binomial_var = sum - sq_sum;

  ExecEnv::log().info( "Genome: {}, Frequency vector size: {}, mean prob: {}, bin mean: {} var: {}, poisson bin mean: {}, var: {}",
                       locus_results.genome, frequency_vector.size(), p, binomial_mean, binomial_var, poisson_binomial_mean, poisson_binomial_var);

  return { frequency_vector, locus_results};

}


