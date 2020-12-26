//
// Created by kellerberrin on 30/11/20.
//

#include "kgl_analysis_mutation_syngen.h"

#include "kel_distribution.h"

#include <sstream>
#include <iomanip>
#include <algorithm>


namespace kgl = kellerberrin::genome;



std::shared_ptr<const kgl::PopulationDB>
kgl::InbreedSynthetic::generateSyntheticPopulation( double lower_inbreeding,
                                                    double upper_inbreeding,
                                                    double step_inbreeding,
                                                    const std::string& super_population,
                                                    const ContigDB& locus_list,
                                                    const LociiVectorArguments& arguments) {
  static std::mutex log_mutex; // logging mutex.

  // Make a synthetic phased diploid population.
  std::shared_ptr<PopulationDB> synthetic_pop_ptr(std::make_shared<PopulationDB>("SyntheticInbreedingPopulation", DataSourceEnum::Genome1000));

  // Entropy source is the Mersenne twister.
  RandomEntropySource entropy_mt;

  // The real unit distribution [1,0] used to draw variants from the locus.
  UniformUnitDistribution unit_distribution;

  // The random boolean used to draw the variant phase (Male or Female).
  RandomBoolean random_boolean;


  // Generate a list of synthetic inbred genomes.
  size_t counter{0};
  std::vector<std::pair<GenomeId_t , double>> inbreeding_vector;
  double inbreeding = lower_inbreeding;
  while(inbreeding <= (upper_inbreeding + 0.000001)) {

    inbreeding_vector.emplace_back(generateSyntheticGenomeId(inbreeding, super_population, counter), inbreeding);

    inbreeding += step_inbreeding;
    ++counter;

  }

  // For all inbred genomes.
  for (auto const& [genome_id, inbreeding_coefficient] : inbreeding_vector) {

    size_t homozygous_count{0};
    size_t heterozygous_count{0};
    std::vector<GenomeId_t> genome_vector {genome_id};

    // Iterate through the locus_list and create each genome with het/hom ratio
    // stochastically defined by the assigned inbreeding coefficient
    for (auto const& [offset, offset_ptr] : locus_list.getMap()) {

      OffsetDBArray variant_vec = offset_ptr->getVariantArray();

      // Generate the minor allele frequencies.
      AlleleFreqVector freq_vector(variant_vec, super_population, arguments.frequencySource());

      // Draw a unit rand and select an allele class.
      double class_selection = unit_distribution.random(entropy_mt.generator());
      // Generate the allele class frequencies.
      AlleleClassFrequencies class_freqs = freq_vector.alleleClassFrequencies(inbreeding_coefficient);
      // Create the synthetic allele mix at the offset.
      switch(freq_vector.selectAlleleClass(class_selection, class_freqs)) {

        case AlleleClassType::MINOR_HOMOZYGOUS: {
          ++homozygous_count;
          // Randomly draw an allele and clone two copies.
          double allele_selection = unit_distribution.random(entropy_mt.generator());
          std::optional<AlleleFreqRecord> selected_allele = freq_vector.selectMinorHomozygous(allele_selection, class_freqs);
          if (selected_allele) {

            std::shared_ptr<Variant> cloned_variant1 = selected_allele->allele()->clone();
            cloned_variant1->updatePhaseId(VariantSequence::DIPLOID_PHASE_A);
            std::shared_ptr<Variant> cloned_variant2 = selected_allele->allele()->clone();
            cloned_variant2->updatePhaseId(VariantSequence::DIPLOID_PHASE_B);
            if (not synthetic_pop_ptr->addVariant( cloned_variant1, genome_vector)) {

              ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                  , genome_id, cloned_variant1->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add female hom variant
            if (not synthetic_pop_ptr->addVariant( cloned_variant2, genome_vector)) {

              ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                  , genome_id, cloned_variant2->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add male hom variant

          } else {

            ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} MINOR_HOMOZYGOUS, could not select minor variant", genome_id);

          }

        }
          break;

        case AlleleClassType::MAJOR_HETEROZYGOUS: {
          ++heterozygous_count;
          // Draw an allele and randomly decide phase.
          double allele_selection = unit_distribution.random(entropy_mt.generator());
          std::optional<AlleleFreqRecord> selected_allele = freq_vector.selectMajorHeterozygous(allele_selection, class_freqs);
          if (selected_allele) {

            std::shared_ptr<Variant> cloned_variant = selected_allele->allele()->clone();
            if (random_boolean.random(entropy_mt.generator())) {

              cloned_variant->updatePhaseId(VariantSequence::DIPLOID_PHASE_A);

            } else {

              cloned_variant->updatePhaseId(VariantSequence::DIPLOID_PHASE_B);

            }
            if (not synthetic_pop_ptr->addVariant( cloned_variant, genome_vector)) {

              ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                  , genome_id, cloned_variant->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          } else {

            ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} MAJOR_HETEROZYGOUS, could not select minor variant", genome_id);

          }

        }
          break;

        case AlleleClassType::MINOR_HETEROZYGOUS: {
          ++heterozygous_count;
          // Draw an allele and clone two copies.
          double allele_selection = unit_distribution.random(entropy_mt.generator());
          std::optional<std::pair<AlleleFreqRecord, AlleleFreqRecord>> selected_alleles = freq_vector.selectMinorHeterozygous(allele_selection, class_freqs);
          if (selected_alleles) {

            std::shared_ptr<Variant> cloned_variant1 = selected_alleles->first.allele()->clone();
            cloned_variant1->updatePhaseId(VariantSequence::DIPLOID_PHASE_A);
            std::shared_ptr<Variant> cloned_variant2 = selected_alleles->second.allele()->clone();
            cloned_variant2->updatePhaseId(VariantSequence::DIPLOID_PHASE_B);
            if (not synthetic_pop_ptr->addVariant( cloned_variant1, genome_vector)) {

              ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                  , genome_id, cloned_variant1->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add female hom variant
            if (not synthetic_pop_ptr->addVariant( cloned_variant2, genome_vector)) {

              ExecEnv::log().error( "InbreedSampling::generateSyntheticPopulation, Genome: {} cannot add variant: {}"
                  , genome_id, cloned_variant2->output(',', VariantOutputIndex::START_0_BASED, false));

            } // add male hom variant

          } else {
            std::scoped_lock log_lock(log_mutex);

            ExecEnv::log().warn( "InbreedSampling::generateSyntheticPopulation, Genome: {} MINOR_HETEROZYGOUS, could not select minor variants", genome_id);
            for (auto const& allele : freq_vector.alleleFrequencies()) {

              ExecEnv::log().warn("InbreedSampling::generateSyntheticPopulation; freq: {}, variant: {}",
                                  allele.frequency(), allele.allele()->output(',', VariantOutputIndex::START_0_BASED, false));

            }

          }

        }
          break;

        case AlleleClassType::MAJOR_HOMOZYGOUS:
          ++homozygous_count;
          // No Action.
          break;

      } // end switch

    } // locii

    ExecEnv::log().info("Synthetic Inbred Genome: {}, Inbreeding: {}, Total: {}, Homozygous: {}, Heterozygous: {}",
                        genome_id, inbreeding_coefficient, (homozygous_count + heterozygous_count), homozygous_count, heterozygous_count);

  } // for all genomes.

  return synthetic_pop_ptr;

}




// Generate an inbreeding encoded synthetic genome
kgl::GenomeId_t kgl::InbreedSynthetic::generateSyntheticGenomeId( double inbreeding,
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
std::pair<bool, double> kgl::InbreedSynthetic::generateInbreeding(const GenomeId_t& genome_id) {

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


