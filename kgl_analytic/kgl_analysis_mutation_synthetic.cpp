//
// Created by kellerberrin on 19/8/20.
//


#include "kel_thread_pool.h"
#include "kel_distribution.h"
#include "kgl_analysis_mutation_inbreed.h"

#include <sstream>
#include <iomanip>


namespace kgl = kellerberrin::genome;

std::shared_ptr<const kgl::DiploidPopulation>
kgl::InbreedingAnalysis::generateSyntheticPopulation( double lower_inbreeding,
                                                      double upper_inbreeding,
                                                      double step_inbreeding,
                                                      const std::string& super_population,
                                                      const std::shared_ptr<const ContigVariant>& locus_list) const {

  std::shared_ptr<DiploidPopulation> synthetic_pop_ptr(std::make_shared<DiploidPopulation>("SyntheticInbreedingPopulation"));

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;
  // The real unit distribution [1,0] used to draw variants from the locus.
  UniformUnitDistribution unit_distribution;
  // The random boolean used to draw the variant phase (Male or Female).
  RandomBoolean random_boolean;
  // The super population AF field
  const std::string AF_field = lookupSuperPopulationField(super_population);

  size_t counter = 0;
  std::vector<std::pair<GenomeId_t , double>> inbreeding_vector;
  double inbreeding = lower_inbreeding;
  while(inbreeding <= (upper_inbreeding + 0.000001)) {

    inbreeding_vector.emplace_back(generateSyntheticGenomeId(inbreeding, super_population, counter), inbreeding);

    inbreeding += step_inbreeding;
    ++counter;

  }

  // Iterate through the locus_list and create each genome with het/hom ratio
  // stochastically defined by the assigned inbreeding coefficient
  for (auto const& [offset, offset_ptr] : locus_list->getMap()) {

    // For each locus get the super-population frequency of each allele
    std::vector<std::pair<std::shared_ptr<const Variant>, double>> variant_freq_vector;
    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      auto [result, AF_value] = processFloatField(variant_ptr, AF_field);
      if (result) {

        variant_freq_vector.emplace_back(variant_ptr, AF_value);

      } else {

        ExecEnv::log().error( "InbreedingAnalysis::generateSyntheticPopulation, Problem reading AF_field: {} for variant: {}"
                             , AF_field, variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

      }

    } // for variants at the offset.

    // For all inbred genomes.
    for (auto const& [genome_id, inbreeding_coefficient] : inbreeding_vector) {

      // Draw a unit rand.
      double random_variant_select = unit_distribution.random(entropy_mt.generator());
      // The offset variant is stochastically selected (or ignored)
      // and then based on the assigned inbreeding coefficient the variant stochastically selected as is heterozygous or homozygous.
      double sum_frequency = 0.0;
      for (auto const& [variant_ptr, frequency] : variant_freq_vector) {

        sum_frequency += frequency;
        if (random_variant_select <= sum_frequency) {
        // Variant has been selected. Draw another unit random to determine if hom or het.
          double random_hom_het = unit_distribution.random(entropy_mt.generator());
          double hom_probability = (inbreeding_coefficient * frequency) + ((1.0 - inbreeding_coefficient) * (frequency * frequency));
          if (random_hom_het <= hom_probability) {
          // Homozygous
            std::vector<GenomeId_t> genome_vector {genome_id};

            std::shared_ptr<Variant> female_variant_ptr = variant_ptr->clone();
            female_variant_ptr->updatePhaseId(VariantSequence::DIPLOID_PHASE_A);
            // Add to the population.
            if (not synthetic_pop_ptr->addVariant( female_variant_ptr, genome_vector)) {

              ExecEnv::log().error( "InbreedingAnalysis::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                                  , genome_id, female_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add female hom variant

            std::shared_ptr<Variant> male_variant_ptr = variant_ptr->clone();
            male_variant_ptr->updatePhaseId(VariantSequence::DIPLOID_PHASE_B);
            // Add to the population.
            if (not synthetic_pop_ptr->addVariant( male_variant_ptr, genome_vector)) {

              ExecEnv::log().error( "InbreedingAnalysis::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                                  , genome_id, male_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add male hom variant

          } else {
          // Heterozygous
          // Randomly select the variant phase and add to the genome.
            std::shared_ptr<Variant> variant_copy_ptr = variant_ptr->clone();
          // Draw a random boolean to determine the variant phase.
            PhaseId_t random_phase = random_boolean.random(entropy_mt.generator()) ? VariantSequence::DIPLOID_PHASE_A : VariantSequence::DIPLOID_PHASE_B;
            variant_copy_ptr->updatePhaseId(random_phase);
          // Add to the population.
            std::vector<GenomeId_t> genome_vector {genome_id};
            if (not synthetic_pop_ptr->addVariant( variant_copy_ptr, genome_vector)) {

              ExecEnv::log().error( "InbreedingAnalysis::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                                  , genome_id, variant_copy_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add het variant

          } // if hom or het

          break; // only select 1 variant.

        } // a variant was selected.

      } // for all offset variants

    } // for all genomes.

  } // for all locus offsets

  for (auto const& [genome_id, genome_ptr] : synthetic_pop_ptr->getMap()) {

    ExecEnv::log().info("Synthetic Inbred Genome: {}, variant count: {}", genome_id, genome_ptr->variantCount());

  }

  return synthetic_pop_ptr;

}

// Generate an inbreeding encoded synthetic genome
kgl::GenomeId_t kgl::InbreedingAnalysis::generateSyntheticGenomeId(double inbreeding, const std::string& super_population, size_t counter) const {

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
std::pair<bool, double> kgl::InbreedingAnalysis::generateInbreeding(const GenomeId_t& genome_id) const {

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