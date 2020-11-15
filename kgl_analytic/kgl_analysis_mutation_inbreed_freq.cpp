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




kgl::AlleleFreqVector::AlleleFreqVector(const std::vector<std::shared_ptr<const Variant>>& variant_vector,
                                        const std::string& frequency_field) {

  // Loop through the variants in the locus..
  for (auto const &variant : variant_vector) {

    auto[result, AF_value] = InbreedSampling::processFloatField(*variant, frequency_field);
    if (not result) {

      // Problem obtaining allele frequency, this allele may not be defined for the specified super population.
      continue;

    } else  {

      // Check for duplicates.
      bool found_duplicate{false};
      for (auto const& allele : allele_frequencies_) {

        if (allele.allele()->analogous(*variant)) {

          found_duplicate = true;
          break;

        }

      }

      if (not found_duplicate) {

        AF_value = std::clamp(AF_value, 0.0, 1.0); // Just in case.
        // Store the variant and super population frequency
        allele_frequencies_.emplace_back(variant, AF_value, frequency_field);

      }

    } // if AF_value

  } // for

}


bool kgl::AlleleFreqVector::checkValidAlleleVector() {

  static std::mutex log_mutex;
  double allele_sum = sumAlleleFrequencies();
  double check_allele_sum = allele_sum - 1.0;
  if (check_allele_sum > epsilon_class_ and check_allele_sum < epsilon_sum_) {
  // If the sum is not too far from 1.0 we adjust the allele frequencies.
  // This handles minor inaccuracies in the allele frequencies in the data files.

    for (auto& allele : allele_frequencies_) {

      double adjusted_frequency = allele.frequency() / allele_sum;
      allele.frequency(adjusted_frequency);

    }

  } else if (check_allele_sum > epsilon_sum_) {
  // The allele frequency error is too big to ignore (arguably).
    std::scoped_lock lock(log_mutex); // Synchronous thread logging.

    std::stringstream ss;
    for (auto const& allele : alleleFrequencies()) {

      ss << allele.frequency() << ", ";

    }
    ExecEnv::log().warn("AlleleFreqVector::checkValidAlleleVector; Sum minor allele freqs: {} > 1.0 ({})",
                        ss.str(), minorAlleleFrequencies());

    for (auto const& allele : alleleFrequencies()) {

      ExecEnv::log().warn("AlleleFreqVector::checkValidAlleleVector; Genome Id: {}, Super Pop: {}, frequency: {}, SNP: {}",
                          allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

    }

    return false;

  }

  return true;

}


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

double kgl::AlleleFreqVector::sumAlleleFrequencies() const {

  double sum_allele_freq{0.0};
  for (auto const& allele_freq : allele_frequencies_) {

    sum_allele_freq += allele_freq.frequency();

  }

  return std::clamp(sum_allele_freq, 0.0, 1.0);

}

double kgl::AlleleFreqVector::minorAlleleFrequencies() const {

  return std::clamp(sumAlleleFrequencies(), 0.0, 1.0);

}

double kgl::AlleleFreqVector::majorAlleleFrequency() const {

  return std::clamp((1.0 - minorAlleleFrequencies()), 0.0, 1.0);

}


double kgl::AlleleFreqVector::majorHomozygous(double inbreeding) const {

  double major_freq = majorAlleleFrequency();
  double prob_hom = (major_freq * inbreeding) + (1.0-inbreeding) * major_freq * major_freq;
  return std::clamp(prob_hom, 0.0, 1.0);

}


// The probability of two minor heterozygous alleles (rare).
double kgl::AlleleFreqVector::minorHeterozygous(double inbreeding) const {

  double sum_hetero_freq{0.0};
  size_t minor_allele_count = allele_frequencies_.size();

  for (size_t idx1 = 0; idx1 < minor_allele_count; ++idx1) {

    for (size_t idx2 = (idx1 + 1); idx2 < minor_allele_count; ++idx2) {

        sum_hetero_freq += 2.0 * allele_frequencies_[idx1].frequency() * allele_frequencies_[idx2].frequency();

    }

  }

  sum_hetero_freq = sum_hetero_freq * (1.0 - inbreeding);
  return std::clamp(sum_hetero_freq, 0.0, 1.0);

}

// The probability of major heterozygous alleles. This is the case where only a single minor allele is recorded (common).
double kgl::AlleleFreqVector::majorHeterozygous(double inbreeding) const {

  double sum_hetero_freq{0.0};
  double major_allele_freq = majorAlleleFrequency();
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_hetero_freq += 2.0 * minor_allele.frequency() * major_allele_freq;

  }

  sum_hetero_freq = sum_hetero_freq * (1.0 - inbreeding);
  return std::clamp(sum_hetero_freq, 0.0, 1.0);

}


// The probability of minor homozygous alleles
double kgl::AlleleFreqVector::minorHomozygous(double inbreeding) const {

  double sum_homozygous_freq{0.0};
  for (auto const& minor_allele  : allele_frequencies_) {

    double minor_freq = minor_allele.frequency();
    sum_homozygous_freq += (minor_freq * inbreeding) + (1.0 - inbreeding) * minor_freq * minor_freq;

  }

  return std::clamp(sum_homozygous_freq, 0.0, 1.0);

}


kgl::AlleleClassFrequencies  kgl::AlleleFreqVector::alleleClassFrequencies(double inbreeding) const {

  static std::mutex log_mutex;

  double major_hom = majorHomozygous(inbreeding);
  double major_het = majorHeterozygous(inbreeding);
  double minor_hom = minorHomozygous(inbreeding);
  double minor_het = minorHeterozygous(inbreeding);

  double sum_freq_classes = major_hom + major_het + minor_hom + minor_het;
  double raw_freq_sum = std::fabs(sum_freq_classes - 1.0);

  if (raw_freq_sum > epsilon_sum_ and inbreeding >= 0.0) {
    std::scoped_lock log_lock(log_mutex);

    ExecEnv::log().error("AlleleFreqVector::alleleClassFrequencies; Invalid Major Hom: {}, Major Het: {}, Minor Hom: {}, Minor Het: {}, Sum: {}, Inbreeding: {}, Allele Count: {}",
                         major_hom, major_het, minor_hom, minor_het, sum_freq_classes, inbreeding, allele_frequencies_.size());

    for (auto const& allele : allele_frequencies_) {

      ExecEnv::log().error("AlleleFreqVector::alleleClassFrequencies; Frequency: {}, Minor Allele: {}",
                           allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

    }

  } else if (raw_freq_sum <= epsilon_sum_) {

    major_hom = major_hom / sum_freq_classes;
    major_het = major_het / sum_freq_classes;
    minor_hom = minor_hom / sum_freq_classes;
    minor_het = minor_het / sum_freq_classes;

  } else if (inbreeding < 0.0 and sum_freq_classes > 0) {

    major_hom = major_hom / sum_freq_classes;
    major_het = major_het / sum_freq_classes;
    minor_hom = minor_hom / sum_freq_classes;
    minor_het = minor_het / sum_freq_classes;

  }

  return AlleleClassFrequencies(major_hom, major_het, minor_hom, minor_het);

}


// Randomly select an allele class outcome based on a unit [0, 1] random number.
kgl::AlleleClassType kgl::AlleleFreqVector::selectAlleleClass(double unit_rand, double inbreeding) const {

  AlleleClassFrequencies class_freqs = alleleClassFrequencies(inbreeding);
  double sum_freqs = class_freqs.minorHeterozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MINOR_HETEROZYGOUS;

  }

  sum_freqs += class_freqs.minorHomozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MINOR_HOMOZYGOUS;

  }

  sum_freqs += class_freqs.majorHomozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MAJOR_HOMOZYGOUS;

  }

  return AlleleClassType::MAJOR_HETEROZYGOUS;

}


std::optional<kgl::AlleleFreqRecord> kgl::AlleleFreqVector::selectMinorAllele(double unit_rand) const {

  if (allele_frequencies_.empty()) {

    ExecEnv::log().error("AlleleFreqVector::selectMinorAllele; attempted to select non-existent minor allele");
    return std::nullopt;

  }

  double freq_sum =  sumAlleleFrequencies();

  if (freq_sum == 0.0) {

    ExecEnv::log().warn("AlleleFreqVector::selectMinorAllele; allele frequencies sum to 0");
    return std::nullopt;

  }

  if (allele_frequencies_.size() == 1) {

    return allele_frequencies_.front();

  }

  for (auto const& allele : allele_frequencies_) {

    if (unit_rand <= (allele.frequency() / freq_sum)) {

      return allele;

    }

  }

  return allele_frequencies_.front();

}


std::optional<std::pair<kgl::AlleleFreqRecord, kgl::AlleleFreqRecord>>
kgl::AlleleFreqVector::selectPairMinorAllele(double unit_rand1, double unit_rand2) const {

  if (allele_frequencies_.size() < 2) {

    ExecEnv::log().error("AlleleFreqVector::selectPairMinorAllele; attempted to select non-existent minor allele pair");
    return std::nullopt;

  }

  double freq_sum =  sumAlleleFrequencies();

  if (freq_sum == 0.0) {

    ExecEnv::log().warn("AlleleFreqVector::selectMinorAllele; allele frequencies sum to 0");
    return std::nullopt;

  }

  if (allele_frequencies_.size() == 2) {

    std::pair<kgl::AlleleFreqRecord, kgl::AlleleFreqRecord> allele_pair{ allele_frequencies_.front(), allele_frequencies_.back() };
    return allele_pair;

  }

  size_t select_index_1{0};
  size_t select_index_2{0};
  double freq_sum_2{0.0};
  std::vector<AlleleFreqRecord> allele_frequencies2;

  size_t index{0};
  bool found_index{false};
  for (auto const& allele : allele_frequencies_) {

    if (unit_rand1 <= (allele.frequency() / freq_sum) and not found_index) {

      select_index_1 = index;
      found_index = true;

    } else {

      ++index;
      freq_sum_2 += allele.frequency();
      allele_frequencies2.push_back(allele);

    }

  }

  if (freq_sum_2 == 0.0) {

    ExecEnv::log().warn("AlleleFreqVector::selectMinorAllele; 2nd allele frequencies sum to 0");
    return std::nullopt;

  }

  index = 0;
  for (auto const& allele : allele_frequencies2) {

    if (unit_rand2 <= (allele.frequency() / freq_sum_2)) {

      select_index_2 = index;
      break;

    }

  }

  std::pair<kgl::AlleleFreqRecord, kgl::AlleleFreqRecord> allele_pair{ allele_frequencies_[select_index_1], allele_frequencies2[select_index_2] };
  return allele_pair;

}


double kgl::AlleleFreqInfo::alleleTypeFrequency(double inbreeding) const {

  switch(alleleType()) {

    case AlleleClassType::MAJOR_HOMOZYGOUS:
      return allele_frequencies_.majorHomozygous(inbreeding);

    case AlleleClassType::MAJOR_HETEROZYGOUS:
      return allele_frequencies_.majorHeterozygous(inbreeding);

    case AlleleClassType::MINOR_HETEROZYGOUS:
      return allele_frequencies_.minorHeterozygous(inbreeding);

    case AlleleClassType::MINOR_HOMOZYGOUS:
      return allele_frequencies_.minorHomozygous(inbreeding);

  }

  return 0.0;

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

    AlleleFreqVector allele_freq_vector(locus_variant_array, gnomad_super_pop);
    if (not allele_freq_vector.checkValidAlleleVector()) {

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
            std::shared_ptr<const Variant> null_variant = allele_freq.allele()->cloneNullVariant();
            AlleleFreqRecord major_allele(null_variant, major_allele_frequency, gnomad_super_pop);
            // Add to the frequency vector.
            frequency_vector.emplace_back(AlleleClassType::MAJOR_HETEROZYGOUS, allele_freq, major_allele, allele_freq_vector);
            break;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset[0]->homozygous(*diploid_offset[1])) {

              frequency_vector.emplace_back(AlleleClassType::MINOR_HOMOZYGOUS, allele_freq, allele_freq, allele_freq_vector);
              break;


            } else {
              // The sample has different alt alleles (rare).
              // Find the second minor allele and obtain it's frequency.
              bool found_second_minor = false;
              size_t second_allele_index{0};
              for (auto const& second_allele_freq : allele_freq_vector.alleleFrequencies()) {

                if (diploid_offset.back()->analogous(*second_allele_freq.allele())) {

                  found_second_minor= true;
                  break;

                }

                ++second_allele_index;

              } // For find second allele frequency.

              if (not found_second_minor) {
                ExecEnv::log().warn("InbreedingCalculation::generateFrequencies; Genome: {}, Not Found Second Minor SNP: {}",
                                    genome_id, diploid_offset[0]->output(',',VariantOutputIndex::START_0_BASED, false));

              } else {

                const AlleleFreqRecord& second_allele = allele_freq_vector.alleleFrequencies().at(second_allele_index);
                frequency_vector.emplace_back(AlleleClassType::MINOR_HETEROZYGOUS, allele_freq, second_allele, allele_freq_vector);
                // No need for further processing. Note that allele_vector is invalid (std::move) at this point.
                break;

              }

            } // Two minor alleles

          }  // if 2 alleles

        } // Valid minor allele frequency.

      } // For all minor alleles in the reference at offset

    } else {
    // No minor alleles.

      if (allele_freq_vector.alleleFrequencies().empty()) {

        ExecEnv::log().error("InbreedingCalculation::generateFrequencies, Genome: {}, zero sized minor allele vector at offset: {}",
                             genome_id, offset);

      } else {

        double major_allele_frequency = allele_freq_vector.majorAlleleFrequency();
        constexpr static const double minimum_major_frequency = 0.01;
        if (major_allele_frequency > minimum_major_frequency) {

          std::shared_ptr<const Variant> null_variant = allele_freq_vector.alleleFrequencies().front().allele()->cloneNullVariant();
          AlleleFreqRecord major_allele(null_variant, major_allele_frequency, gnomad_super_pop);
          frequency_vector.emplace_back(AlleleClassType::MAJOR_HOMOZYGOUS, major_allele, major_allele, allele_freq_vector);

        }

      }

    }

  } // For all offset locii.

  // Generate some frequency statistics.
  locus_results.total_allele_count = frequency_vector.size();
  for (auto const& allele_freq : frequency_vector) {

    // Get the allele class frequencies.
    AlleleClassFrequencies class_frequencies = allele_freq.alleleFrequencies().alleleClassFrequencies(0.0);
    locus_results.major_homo_freq += class_frequencies.majorHomozygous();
    locus_results.minor_homo_freq += class_frequencies.minorHomozygous();
    locus_results.major_hetero_freq += class_frequencies.majorHeterozygous();
    locus_results.minor_hetero_freq += class_frequencies.minorHeterozygous();

    // Count the actual allele classes.
    switch(allele_freq.alleleType()) {

      case AlleleClassType::MINOR_HOMOZYGOUS:
        ++locus_results.minor_homo_count;
        break;

      case AlleleClassType::MAJOR_HETEROZYGOUS:
        ++locus_results.major_hetero_count;
        break;

      case AlleleClassType::MINOR_HETEROZYGOUS:
        ++locus_results.minor_hetero_count;
        break;

      case AlleleClassType::MAJOR_HOMOZYGOUS:
        ++locus_results.major_homo_count;
        break;

    }

  }

  return { frequency_vector, locus_results };

}


