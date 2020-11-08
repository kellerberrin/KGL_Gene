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


bool kgl::AlleleFreqVector::checkDuplicates() const {

  size_t minor_allele_count = allele_frequencies_.size();

  for (size_t idx1 = 0; idx1 < minor_allele_count; ++idx1) {

    for (size_t idx2 = (idx1 + 1); idx2 < minor_allele_count; ++idx2) {

      if (allele_frequencies_[idx1].allele()->analogous(*allele_frequencies_[idx1].allele())) {

        return true;

      }

    }

  }

  return false;

}


double kgl::AlleleFreqVector::majorAlleleFrequency() const {

  double sum_allele_freq{0.0};
  for (auto const& allele_freq : allele_frequencies_) {

    sum_allele_freq += allele_freq.frequency();

  }

  return std::clamp((1.0 - sum_allele_freq), 0.0, 1.0);

}


double kgl::AlleleFreqVector::majorHomozygous() const {


  return majorAlleleFrequency() * majorAlleleFrequency();

}


// The probability of two minor heterozygous alleles (rare).
double kgl::AlleleFreqVector::minorHeterozygous() const {

  double sum_hetero_freq{0.0};
  size_t minor_allele_count = allele_frequencies_.size();

  for (size_t idx1 = 0; idx1 < minor_allele_count; ++idx1) {

    for (size_t idx2 = (idx1 + 1); idx2 < minor_allele_count; ++idx2) {

        sum_hetero_freq += 2.0 * allele_frequencies_[idx1].frequency() * allele_frequencies_[idx2].frequency();

    }

  }

  return sum_hetero_freq;

}

// The probability of major heterozygous alleles. This is the case where only a single minor allele is recorded (common).
double kgl::AlleleFreqVector::majorHeterozygous() const {

  double sum_hetero_freq{0.0};
  double major_allele_freq = majorAlleleFrequency();
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_hetero_freq += 2.0 * minor_allele.frequency() * major_allele_freq;

  }

  return sum_hetero_freq;

}


// The probability of minor homozygous alleles
double kgl::AlleleFreqVector::minorHomozygous() const {

  double sum_homozygous_freq{0.0};
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_homozygous_freq += minor_allele.frequency() * minor_allele.frequency();

  }

  return sum_homozygous_freq;

}



double kgl::AlleleFreqVector::sumFrequencies() const {

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate Allele Frequencies using Gnomad.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<kgl::AlleleFreqInfo>, kgl::LocusResults>
kgl::InbreedingCalculation::generateFrequencies(const GenomeId_t& genome_id,
                                                const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                const std::string& super_population_field,
                                                const std::shared_ptr<const ContigVariant>& locus_list) {

  static std::mutex log_mutex;  // Synchronous thread logging.
  std::vector<AlleleFreqInfo> frequency_vector;
  LocusResults locus_results;
  locus_results.genome = genome_id;

  // Convert to the super population field codes used in Gnomad.
  const std::string gnomad_super_pop = InbreedSampling::lookupSuperPopulationField(super_population_field);
  // Diploid contig, only want SNP variants.
  auto snp_contig_ptr = contig_ptr->filterVariants(SNPFilter());

  // For all offsets.
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // Get the allele frequencies.
   auto locus_variant_array = offset_ptr->getVariantArray();

    std::vector<AlleleFreqRecord> allele_vector;
    // Loop through the variants in the locus..
    for (auto const &locus_variant : locus_variant_array) {

      auto[result, AF_value] = InbreedSampling::processFloatField(*locus_variant, gnomad_super_pop);
      if (not result) {

        // Problem obtaining allele frequency.
        ExecEnv::log().warn("InbreedingAnalysis::generateFrequencies; Genome: {}, Unable to obtain Locus allele frequency: {} for SNP: {}",
                            genome_id, super_population_field, locus_variant->output(',', VariantOutputIndex::START_0_BASED, false));

      } else  {

        AF_value = std::clamp(AF_value, 0.0, 1.0); // Just in case.
        // Store the Gnomad variant and super population frequency
        allele_vector.emplace_back(locus_variant, AF_value);

      }

    }

    AlleleFreqVector allele_freq_vector(allele_vector);

    if (allele_freq_vector.checkDuplicates()) {
      std::scoped_lock lock(log_mutex); // Synchronous thread logging.

      for (auto const &allele : allele_freq_vector.alleleFrequencies()) {

        ExecEnv::log().warn("InbreedingCalculation::generateFrequencies; Genome Id: {}, Super Pop: {}, frequency: {}, Duplicate SNP: {}",
                            genome_id, super_population_field, allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

      }

      continue; // next locus.

    }

    if (allele_freq_vector.sumFrequencies() > 1) {
      std::scoped_lock lock(log_mutex); // Synchronous thread logging.

      std::stringstream ss;
      for (auto const& allele : allele_freq_vector.alleleFrequencies()) {

        ss << allele.frequency() << ", ";

      }
      ExecEnv::log().warn("InbreedingCalculation::generateFrequencies; Genome Id: {}, Super Pop: {}, Sum minor allele freqs: {} > 1.0 ({})",
                          genome_id, super_population_field, ss.str(), allele_freq_vector.sumFrequencies());

      for (auto const& allele : allele_freq_vector.alleleFrequencies()) {

        ExecEnv::log().warn("InbreedingCalculation::generateFrequencies; Genome Id: {}, Super Pop: {}, frequency: {}, SNP: {}",
                            genome_id, super_population_field, allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

      }

      continue; // Next locus.

    }


    // Check the diploid genome for any minor alleles at this location.
    auto diploid_variant_opt = snp_contig_ptr->findOffsetArray(offset);
    // If minor alleles are at the diploid location.
    if (diploid_variant_opt) {

      auto const &diploid_offset = diploid_variant_opt.value();

      // Loop through the reference frequency vector
      for (auto const &allele_freq : allele_freq_vector.alleleFrequencies()) {
        // The diploid offset can have 1 or 2 minor alleles
        // We need only check the first against the locus which is assumed to list all possible minor alleles.
        if (diploid_offset.front()->analogous(*allele_freq.allele())) {

          if (diploid_offset.size() == 1) {
            // The sample is minor allele heterozygous
            // Calculate the major allele frequency as the complement of the sum of minor alleles.
            double major_allele_frequency = allele_freq_vector.majorAlleleFrequency();
            // Add to the frequency vector.
            frequency_vector.emplace_back(MinorAlleleType::MAJOR_HETEROZYGOUS, allele_freq.frequency(), major_allele_frequency, allele_freq_vector);
            break;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset[0]->homozygous(*diploid_offset[1])) {

              frequency_vector.emplace_back(MinorAlleleType::MINOR_HOMOZYGOUS, allele_freq.frequency(), allele_freq.frequency(), allele_freq_vector);
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

                frequency_vector.emplace_back(MinorAlleleType::MINOR_HETEROZYGOUS, allele_freq.frequency(), second_frequency, allele_freq_vector);
                // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
                break;

              }

            } // Two minor alleles

          }  // if 2 alleles

        } // Valid minor allele frequency.

      } // For all minor alleles in the reference at offset

    } else {
    // No minor alleles.

      double major_allele_frequency = allele_freq_vector.majorAlleleFrequency();
      frequency_vector.emplace_back(MinorAlleleType::MAJOR_HOMOZYGOUS, major_allele_frequency, major_allele_frequency, allele_freq_vector);

    }

  } // For all offset locii.

  // Generate some frequency statistics.
  locus_results.total_allele_count = frequency_vector.size();
  for (auto const& allele_freq : frequency_vector) {

    locus_results.major_homo_freq += allele_freq.alleleFrequencies().majorHomozygous();
    locus_results.minor_homo_freq += allele_freq.alleleFrequencies().minorHomozygous();
    locus_results.major_hetero_freq += allele_freq.alleleFrequencies().majorHeterozygous();
    locus_results.minor_hetero_freq += allele_freq.alleleFrequencies().minorHeterozygous();

    switch(allele_freq.alleleType()) {

      case MinorAlleleType::MINOR_HOMOZYGOUS:
        ++locus_results.minor_homo_count;
        break;

      case MinorAlleleType::MAJOR_HETEROZYGOUS:
        ++locus_results.major_hetero_count;
        break;

      case MinorAlleleType::MINOR_HETEROZYGOUS:
        ++locus_results.minor_hetero_count;
        break;

      case MinorAlleleType::MAJOR_HOMOZYGOUS:
        ++locus_results.major_homo_count;
        break;

    }

  }

  return { frequency_vector, locus_results};

}


