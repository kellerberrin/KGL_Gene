//
// Created by kellerberrin on 19/8/20.
//


#include "kel_distribution.h"
#include "kgl_analysis_mutation_inbreed_aux.h"

#include <sstream>
#include <iomanip>
#include <algorithm>


namespace kgl = kellerberrin::genome;




std::shared_ptr<const kgl::DiploidPopulation>
kgl::InbreedSampling::generateSyntheticPopulation( double lower_inbreeding,
                                                   double upper_inbreeding,
                                                   double step_inbreeding,
                                                   const std::string& super_population,
                                                   const ContigVariant& locus_list) {

  std::shared_ptr<DiploidPopulation> synthetic_pop_ptr(std::make_shared<DiploidPopulation>("SyntheticInbreedingPopulation"));

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;
  // The real unit distribution [1,0] used to draw variants from the locus.
  UniformUnitDistribution unit_distribution;

  // The random boolean used to draw the variant phase (Male or Female).
  RandomBoolean random_boolean;
  // The super population AF field
  const std::string AF_field = InbreedSampling::lookupSuperPopulationField(super_population);

  // Generate a list of synthetic inbred genomes.
  size_t counter = 0;
  std::vector<std::pair<GenomeId_t , double>> inbreeding_vector;
  double inbreeding = lower_inbreeding;
  while(inbreeding <= (upper_inbreeding + 0.000001)) {

    inbreeding_vector.emplace_back(generateSyntheticGenomeId(inbreeding, super_population, counter), inbreeding);

    inbreeding += step_inbreeding;
    ++counter;

  }

  // For all inbred genomes.
  for (auto const& [genome_id, inbreeding_coefficient] : inbreeding_vector) {

    size_t homozygous_count = 0;
    size_t heterozygous_count = 0;

  // Iterate through the locus_list and create each genome with het/hom ratio
  // stochastically defined by the assigned inbreeding coefficient
    for (auto const& [offset, offset_ptr] : locus_list.getMap()) {

      auto variant_vec = offset_ptr->getVariantArray();
      // Draw a unit rand.
      double random_variant_selection = unit_distribution.random(entropy_mt.generator());
      double probability_sum = 0;
      std::pair<std::shared_ptr<const Variant>, double> variant_freq{nullptr, 0.0};
      std::vector<std::pair<std::shared_ptr<const Variant>, double>> minor_allele_vec;
      bool minor_allele_selected = false;
      for (auto const& variant : variant_vec) {

        auto [result, AF_value] = InbreedSampling::processFloatField(*variant, AF_field);
        if (result) {

          probability_sum += AF_value;
          if (random_variant_selection <= probability_sum and not minor_allele_selected) {

            variant_freq.first = variant;
            variant_freq.second = AF_value;
            minor_allele_selected = true;

          } else {

            minor_allele_vec.emplace_back(variant, AF_value);

          }

        }

      }

// If we don't have a variant selected then continue
      if (not minor_allele_selected) {

        continue;

      }

      // Given the selection of the minor allele, calc the conditional probability of a homozygous allele
      double hom_conditional_probability = inbreeding_coefficient + ((1.0 - inbreeding_coefficient) * variant_freq.second);
      hom_conditional_probability = std::clamp<double>(hom_conditional_probability, 0.0, 1.0); // limit to [0, 1]
      // Draw for the selection.
      double random_hom_het = unit_distribution.random(entropy_mt.generator());
      if (random_hom_het <= hom_conditional_probability) {
        // Homozygous
        ++homozygous_count;

        std::vector<GenomeId_t> genome_vector {genome_id};
        std::shared_ptr<Variant> female_variant_ptr = variant_freq.first->clone();
        female_variant_ptr->updatePhaseId(VariantSequence::DIPLOID_PHASE_A);
        // Add to the population.
        if (not synthetic_pop_ptr->addVariant( female_variant_ptr, genome_vector)) {

          ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
              , genome_id, female_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        } // add female hom variant

        std::shared_ptr<Variant> male_variant_ptr = variant_freq.first->clone();
        male_variant_ptr->updatePhaseId(VariantSequence::DIPLOID_PHASE_B);
        // Add to the population.
        if (not synthetic_pop_ptr->addVariant( male_variant_ptr, genome_vector)) {

          ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
              , genome_id, male_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        } // add male hom variant


      } else {
        // Heterozygous
        ++heterozygous_count;
        // Randomly select the variant phase and add to the genome.
        std::shared_ptr<Variant> variant_copy_ptr = variant_freq.first->clone();
        // Draw a random boolean to determine the variant phase.
        PhaseId_t allele_phase;
        PhaseId_t alt_allele_phase;
        if (random_boolean.random(entropy_mt.generator())) {

          allele_phase = VariantSequence::DIPLOID_PHASE_A;
          alt_allele_phase = VariantSequence::DIPLOID_PHASE_B;

        } else {

          alt_allele_phase = VariantSequence::DIPLOID_PHASE_A;
          allele_phase = VariantSequence::DIPLOID_PHASE_B;

        }

        variant_copy_ptr->updatePhaseId(allele_phase);
        // Add to the population.
        std::vector<GenomeId_t> genome_vector {genome_id};
        if (not synthetic_pop_ptr->addVariant( variant_copy_ptr, genome_vector)) {

          ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
              , genome_id, variant_copy_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        } // add het variant

        // Check if the other allele is one of the minor alleles (if they exist).
        if (not minor_allele_vec.empty() and hom_conditional_probability < 1.0) {

          double sum_alt_allele_prob = 0.0;
          double alt_allele_select = unit_distribution.random(entropy_mt.generator());
          for (auto const& [minor_allele, frequency] : minor_allele_vec) {

            double het_conditional_probability = 2.0 * (1.0 - inbreeding_coefficient) * frequency;
            sum_alt_allele_prob += het_conditional_probability;
            if (alt_allele_select <= sum_alt_allele_prob) {

              std::shared_ptr<Variant> cloned_allele_ptr = minor_allele->clone();
              cloned_allele_ptr->updatePhaseId(alt_allele_phase);
              // Add to the population.
              if (not synthetic_pop_ptr->addVariant( cloned_allele_ptr, genome_vector)) {

                ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                    , genome_id, cloned_allele_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

              } // add minor allele variant

              break;

            }

          }

        }

      } // if hom or het

    } // locii

    ExecEnv::log().info("Synthetic Inbred Genome: {}, Inbreeding: {}, Total: {}, Homozygous: {}, Heterozygous: {}",
                        genome_id, inbreeding_coefficient, (homozygous_count + heterozygous_count), homozygous_count, heterozygous_count);

  } // for all genomes.

  return synthetic_pop_ptr;

}

// Generate an inbreeding encoded synthetic genome
kgl::GenomeId_t kgl::InbreedSampling::generateSyntheticGenomeId( double inbreeding,
                                                                 const std::string& super_population,
                                                                 size_t counter) {

  // Generate inbreeding vector and create genome ids.
  std::stringstream genome_id_stream;
  genome_id_stream << std::fixed;
  genome_id_stream << std::setprecision(0);

  if (inbreeding >= 0) {

    genome_id_stream << super_population << "_" << (inbreeding * SYNTHETIC_GENOME) << "_" << counter;

  } else {

    genome_id_stream << super_population << "_N" << (-inbreeding * SYNTHETIC_GENOME) << "_" << counter;

  }

  return genome_id_stream.str();

}

// Recreate the inbreeding coefficient from the synthetic genome id.
std::pair<bool, double> kgl::InbreedSampling::generateInbreeding(const GenomeId_t& genome_id) {

  // Decode the inbreeding coefficient.
  bool valid_value = false;
  bool negative = false;
  double inbreed_coefficient = -1000.0; // Impossible value.
  auto first_pos = genome_id.find_first_of("N");
  if (first_pos == std::string::npos) {

    first_pos = genome_id.find_first_of("_");

  } else {

    negative = true;

  }

  if (first_pos != std::string::npos) {

    ++first_pos;
    auto second_pos = genome_id.find_first_of("_", first_pos);
    if (second_pos != std::string::npos) {

      std::string coefficient_string = genome_id.substr( first_pos, second_pos - first_pos);
      try {

        size_t coefficent = std::stod(coefficient_string);
        inbreed_coefficient = static_cast<double>(coefficent)/static_cast<double>(SYNTHETIC_GENOME);
        valid_value = true;
        if (negative) {

          inbreed_coefficient = -1.0 * inbreed_coefficient;

        }

      }
      catch(std::exception& e) {

        ExecEnv::log().warn("InbreedingAnalysis::generateInbreeding, Genome id: {}, cannot convert coefficient string: {}, reason: {}",
                            genome_id, coefficient_string, e.what());

      }

    }

  }

  return {valid_value, inbreed_coefficient};

}