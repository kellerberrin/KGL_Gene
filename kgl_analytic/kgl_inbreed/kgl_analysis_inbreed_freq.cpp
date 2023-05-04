//
// Created by kellerberrin on 2/11/20.
//

#include "kgl_analysis_inbreed_locus.h"
#include "kgl_variant_filter.h"
#include "kgl_analysis_inbreed_calc.h"
#include "kel_distribution.h"

#include <fstream>
#include <algorithm>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;



kgl::AlleleFreqVector::AlleleFreqVector(const OffsetDBArray& variant_vector,
                                        const std::string& frequency_field) {

  // Loop through the variants in the locus..
  for (auto const &variant : variant_vector) {

    auto opt_value = FrequencyDatabaseRead::superPopFrequency(*variant, frequency_field);
    if (not opt_value) {

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

        double variant_freq = std::clamp(opt_value.value(), 0.0, 1.0); // Just in case.
        // Store the variant and super population frequency
        allele_frequencies_.emplace_back(variant, variant_freq);

      }

    } // if AF_value

  } // for

}



bool kgl::AlleleFreqVector::checkValidAlleleVector() {

  double allele_sum = sumAlleleFrequencies();
  double check_allele_sum = allele_sum - 1.0;
  // If greater than 1, then reject.
  if (check_allele_sum > epsilon_class_) {

    return false;

  }

  // Minor allele array must be non-empty
  return not allele_frequencies_.empty();

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

  return sum_allele_freq;

}

double kgl::AlleleFreqVector::minorAlleleFrequencies() const {

  return std::clamp(sumAlleleFrequencies(), 0.0, 1.0);

}

double kgl::AlleleFreqVector::majorAlleleFrequency() const {

  return std::clamp((1.0 - minorAlleleFrequencies()), 0.0, 1.0);

}



kgl::AlleleClassFrequencies  kgl::AlleleFreqVector::unadjustedAlleleClassFrequencies(double inbreeding) const {

  static std::mutex log_mutex;
  std::vector<double> minor_allele_frequencies;

  double sum_minor_freq{0.0};
  for (auto const& minor_allele  : allele_frequencies_) {

    sum_minor_freq += minor_allele.frequency();

  }

  // Can be data problems with allele frequencies summing > 1.0
  double major_frequency = std::max(0.0, (1.0 - sum_minor_freq));
  for (auto const& minor_allele  : allele_frequencies_) {

    if (sum_minor_freq > 1.0) {

      minor_allele_frequencies.push_back((minor_allele.frequency() / sum_minor_freq));

    } else {

      minor_allele_frequencies.push_back(minor_allele.frequency());

    }

  }

  double minor_homozygous{0.0};
  for (auto const& minor_frequency  : minor_allele_frequencies) {

    minor_homozygous += (inbreeding * minor_frequency) + ((1.0 - inbreeding) * minor_frequency * minor_frequency);

  }


  double minor_heterozygous{0.0};
  size_t minor_allele_count = minor_allele_frequencies.size();

  for (size_t idx1 = 0; idx1 < minor_allele_count; ++idx1) {

    for (size_t idx2 = (idx1+1); idx2 < minor_allele_count; ++idx2) {

      minor_heterozygous += (1.0 - inbreeding) * 2.0 * minor_allele_frequencies[idx1] * minor_allele_frequencies[idx2];

    }

  }

  double major_homozygous = (inbreeding * major_frequency) + ((1.0 - inbreeding) * major_frequency * major_frequency);

  double major_heterozygous{0.0};
  for (auto const& minor_frequency  : minor_allele_frequencies) {

    major_heterozygous += (1.0 - inbreeding) * 2.0 * major_frequency * minor_frequency;

  }


  double sum_freq_classes = major_homozygous + major_heterozygous + minor_homozygous + minor_heterozygous;
  double raw_freq_sum = std::fabs(sum_freq_classes - 1.0);

  if (raw_freq_sum > epsilon_class_) {
    std::scoped_lock log_lock(log_mutex);

    ExecEnv::log().error("AlleleFreqVector::alleleClassFrequencies; Invalid Major Hom: {}, Major Het: {}, Minor Hom: {}, Minor Het: {}, Sum: {}, Inbreeding: {}, Allele Count: {}",
                         major_homozygous, major_heterozygous, minor_homozygous, minor_heterozygous, sum_freq_classes, inbreeding, allele_frequencies_.size());

    for (auto const& allele : allele_frequencies_) {

      ExecEnv::log().error("AlleleFreqVector::alleleClassFrequencies; Frequency: {}, Minor Allele: {}",
                           allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

    }

  }

  return AlleleClassFrequencies(major_homozygous, major_heterozygous, minor_homozygous, minor_heterozygous, inbreeding);

}


kgl::AlleleClassFrequencies  kgl::AlleleFreqVector::alleleClassFrequencies(double inbreeding) const {


  AlleleClassFrequencies  class_freqs = unadjustedAlleleClassFrequencies(inbreeding);

  class_freqs.normalize();

  return class_freqs;

}


// Randomly select an allele class outcome based on a unit [0, 1] random number.
kgl::AlleleClassType kgl::AlleleFreqVector::selectAlleleClass(double unit_rand, const AlleleClassFrequencies& class_freqs) const {

  double sum_freqs = class_freqs.minorHomozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MINOR_HOMOZYGOUS;

  }

  sum_freqs += class_freqs.minorHeterozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MINOR_HETEROZYGOUS;

  }

  sum_freqs += class_freqs.majorHomozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MAJOR_HOMOZYGOUS;

  }

  sum_freqs += class_freqs.majorHeterozygous();

  if (unit_rand <= sum_freqs) {

    return AlleleClassType::MAJOR_HETEROZYGOUS;

  }

  ExecEnv::log().warn("AlleleFreqVector::selectAlleleClass; Class select error, Inbreeding, MajorHom: {}, MajorHet: {}, MinorHom: {}, MinorHet {}",
                      class_freqs.inbreeding(), class_freqs.majorHomozygous(), class_freqs.majorHeterozygous(),
                      class_freqs.minorHomozygous(), class_freqs.minorHeterozygous());

  return AlleleClassType::MAJOR_HOMOZYGOUS;

}


std::optional<kgl::AlleleFreqRecord> kgl::AlleleFreqVector::selectMinorHomozygous(double unit_rand,
                                                                                  const AlleleClassFrequencies& class_freqs) const {

  if (allele_frequencies_.empty()) {

    ExecEnv::log().error("AlleleFreqVector::selectMinorHomozygous; attempted to select non-existent minor allele");
    return std::nullopt;

  }

  if (class_freqs.minorHomozygous() == 0.0) {

    ExecEnv::log().warn("AlleleFreqVector::selectMinorHomozygous; minor homozygous frequency is 0");
    return std::nullopt;

  }

  if (allele_frequencies_.size() == 1) {

    return allele_frequencies_.front();

  }

  double allele_freq_sum{0.0};
  for (auto const& allele : allele_frequencies_) {

    double allele_freq = allele.frequency();
    double hom_prob = (allele_freq * class_freqs.inbreeding()) + (1.0 - class_freqs.inbreeding()) * allele_freq * allele_freq;
    allele_freq_sum += hom_prob / class_freqs.minorHomozygous();
    if (unit_rand <= allele_freq_sum) {

      return allele;

    }

  }

  ExecEnv::log().warn("AlleleFreqVector::selectMinorHomozygous; did not select allele");
  return std::nullopt;

}


// Randomly select a major heterozygous allele based on a unit [0,1] random number, std::nullopt if error (no minor allele).
std::optional<kgl::AlleleFreqRecord> kgl::AlleleFreqVector::selectMajorHeterozygous(double unit_rand,
                                                                                    const AlleleClassFrequencies& class_freqs) const {

  if (allele_frequencies_.empty()) {

    ExecEnv::log().error("AlleleFreqVector::selectMajorHeterozygous; attempted to select non-existent minor allele");
    return std::nullopt;

  }

  if (class_freqs.majorHeterozygous() == 0.0) {

    ExecEnv::log().warn("AlleleFreqVector::selectMajorHeterozygous; minor homozygous frequency is 0");
    return std::nullopt;

  }

  if (allele_frequencies_.size() == 1) {

    return allele_frequencies_.front();

  }

  double allele_freq_sum{0.0};
  double major_freq = majorAlleleFrequency();
  for (auto const& allele : allele_frequencies_) {

    double allele_freq = allele.frequency();
    double het_prob = (1.0 - class_freqs.inbreeding()) * 2.0 * major_freq * allele_freq;
    allele_freq_sum += het_prob / class_freqs.majorHeterozygous();
    if (unit_rand <= allele_freq_sum) {

      return allele;

    }

  }

  static std::mutex log_mutex;
  std::scoped_lock log_lock(log_mutex);

  ExecEnv::log().warn("AlleleFreqVector::selectMajorHeterozygous; did not select allele, rand: {}, sum: {}, major het freq: {}");
  for (auto const& allele : allele_frequencies_) {

    ExecEnv::log().warn("AlleleFreqVector::selectMajorHeterozygous; allele frequency {}, allele: {}",
                        allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

  }
  return std::nullopt;

}


// Randomly select a pair of distinct minor alleles based on two random numbers, std::nullopt if error (not two minor alleles).
std::optional<std::pair<kgl::AlleleFreqRecord, kgl::AlleleFreqRecord>>
kgl::AlleleFreqVector::selectMinorHeterozygous(double unit_rand, const AlleleClassFrequencies& class_freqs) const {

  if (allele_frequencies_.size() < 2) {

    ExecEnv::log().error("AlleleFreqVector::selectMinorHeterozygous; attempted to select non-existent minor alleles");
    return std::nullopt;

  }

  if (class_freqs.minorHeterozygous() == 0.0) {

    ExecEnv::log().warn("AlleleFreqVector::selectMinorHeterozygous; minor homozygous frequency is 0");
    return std::nullopt;

  }

  if (allele_frequencies_.size() == 2) {

    std::pair<kgl::AlleleFreqRecord, kgl::AlleleFreqRecord> pair_alleles{ allele_frequencies_.front(), allele_frequencies_.back() };
    return pair_alleles;

  }

  double allele_freq_sum{0.0};
  size_t allele_count = allele_frequencies_.size();
  for (size_t idx1 = 0; idx1 < allele_count; ++idx1) {

    double allele1_freq = allele_frequencies_[idx1].frequency();

    for (size_t idx2 = (idx1 + 1); idx2 < allele_count; ++idx2) {

      double allele2_freq = allele_frequencies_[idx2].frequency();
      double het_prob = (1.0 - class_freqs.inbreeding()) * 2.0 * allele1_freq * allele2_freq;
      allele_freq_sum += het_prob / class_freqs.minorHeterozygous();
      if (unit_rand <= allele_freq_sum) {

        std::pair<kgl::AlleleFreqRecord, kgl::AlleleFreqRecord> pair_alleles{ allele_frequencies_[idx1], allele_frequencies_[idx2] };
        return pair_alleles;

      }

    }

  }

  ExecEnv::log().warn("AlleleFreqVector::selectMinorHeterozygous; did not select alleles");
  return std::nullopt;

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
                                                const std::shared_ptr<const ContigDB>& contig_ptr,
                                                const std::string& super_population_field,
                                                const std::shared_ptr<const ContigDB>& locus_list) {

  std::vector<AlleleFreqInfo> frequency_vector;
  LocusResults locus_results;
  locus_results.genome = genome_id;

  // Diploid contig, only want SNP variants.
  auto snp_contig_ptr = contig_ptr->viewFilter(SNPFilter());

  // For all offsets.
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // Get the allele frequencies.
    const OffsetDBArray& locus_variant_array = offset_ptr->getVariantArray();

    AlleleFreqVector allele_freq_vector(locus_variant_array, super_population_field);
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
            AlleleFreqRecord major_allele(null_variant, major_allele_frequency);
            // Add to the frequency vector.
            frequency_vector.emplace_back(AlleleClassType::MAJOR_HETEROZYGOUS, allele_freq, major_allele, allele_freq_vector);
            break;

          } else if (diploid_offset.size() == 2) {

            if (diploid_offset.front()->homozygous(*diploid_offset.back())) {

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
                                    genome_id, diploid_offset.front()->output(',',VariantOutputIndex::START_0_BASED, false));

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
          AlleleFreqRecord major_allele(null_variant, major_allele_frequency);
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


