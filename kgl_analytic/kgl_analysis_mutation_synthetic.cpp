//
// Created by kellerberrin on 19/8/20.
//


#include "kel_distribution.h"
#include "kgl_analysis_mutation_inbreed_aux.h"

#include <sstream>
#include <iomanip>


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
    double expected_homozygous = 0.0;
    double expected_variant = 0.0;

  // Iterate through the locus_list and create each genome with het/hom ratio
  // stochastically defined by the assigned inbreeding coefficient
    for (auto const& [offset, offset_ptr] : locus_list.getMap()) {

      // For each locus select a variant based on the allele frequency.
      std::pair<std::shared_ptr<const Variant>, double> variant_freq{nullptr, 0.0};
      auto variant_vec = offset_ptr->getVariantArray();
      if (variant_vec.empty()) {

        continue;

      } else if (variant_vec.size() == 1) {

        variant_freq.first = variant_vec.front();

      } else {

        UniformIntegerDistribution variant_index(0, variant_vec.size()-1);
        size_t drawn_index = variant_index.random(entropy_mt.generator());
        variant_freq.first = variant_vec[drawn_index];

      }
      bool variant_selected = false;
      double random_variant_selection = unit_distribution.random(entropy_mt.generator());
      auto [result, AF_value] = InbreedSampling::processFloatField(*variant_freq.first, AF_field);
      if (result) {

        expected_variant += AF_value;
        if (AF_value >= random_variant_selection) {

          variant_selected = true;
          variant_freq.second = AF_value;
          if (AF_value == 0.0) {

            ExecEnv::log().warn( "InbreedSampling::generateSyntheticPopulation, Zero AF_field: {}, random: {}, variant: {}"
                , AF_field, random_variant_selection, variant_vec[0]->output(',', VariantOutputIndex::START_0_BASED, false));

          }

        }

      } else {

        ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Problem reading AF_field: {} for variant: {}"
                             , AF_field, variant_vec[0]->output(',', VariantOutputIndex::START_0_BASED, false));
        continue;

      }

      auto const& [variant_ptr, frequency] = variant_freq;
      if (not (variant_ptr and variant_selected)) {

        continue;

      }

      double random_hom_het = unit_distribution.random(entropy_mt.generator());
      double hom_probability = (inbreeding_coefficient * frequency) + ((1.0 - inbreeding_coefficient) * (frequency * frequency));
      expected_homozygous += hom_probability >= 0.0 ? hom_probability : 0.0;
      if (random_hom_het <= hom_probability) {
        // Homozygous
        std::vector<GenomeId_t> genome_vector {genome_id};

        std::shared_ptr<Variant> female_variant_ptr = variant_ptr->clone();
        female_variant_ptr->updatePhaseId(VariantSequence::DIPLOID_PHASE_A);
          // Add to the population.
        if (not synthetic_pop_ptr->addVariant( female_variant_ptr, genome_vector)) {

          ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                              , genome_id, female_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        } // add female hom variant

        std::shared_ptr<Variant> male_variant_ptr = variant_ptr->clone();
        male_variant_ptr->updatePhaseId(VariantSequence::DIPLOID_PHASE_B);
          // Add to the population.
        if (not synthetic_pop_ptr->addVariant( male_variant_ptr, genome_vector)) {

          ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                              , genome_id, male_variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        } // add male hom variant
        ++homozygous_count;


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

          ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                              , genome_id, variant_copy_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

        } // add het variant

        ++heterozygous_count;

      } // if hom or het

    } // locii

    ExecEnv::log().info("Synthetic Inbred Genome: {}, Inbreeding: {}, Total: {}, Expected Total: {}, Homozygous: {}, ExpectedHom: {}, Heterozygous: {}",
                        genome_id, inbreeding_coefficient, (homozygous_count + heterozygous_count), expected_variant, homozygous_count, expected_homozygous, heterozygous_count);

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