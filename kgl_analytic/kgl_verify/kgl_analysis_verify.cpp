//
// Created by kellerberrin on 7/11/20.
//

#include "kgl_analysis_verify.h"
#include "kgl_variant_filter.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::VerifyAnalysis::initializeAnalysis( const std::string& work_directory,
                                              const ActiveParameterList& named_parameters,
                                              std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::VerifyAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  ExecEnv::log().info("VCF File Read for Analysis Id: {} called with Variant Population", ident());

  auto file_characteristic = data_ptr->dataCharacteristic();

  // Check if we have a population.
  if (file_characteristic.data_implementation == DataImplEnum::PopulationVariant) {

    std::shared_ptr<const PopulationDB> population = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

    // Apply some filters inSitu, note this modifies the population structure so we need to cast away const-ness
    auto non_const_population = std::const_pointer_cast<PopulationDB>(population);

    if (file_characteristic.data_source == DataSourceEnum::Genome1000
        or file_characteristic.data_source == DataSourceEnum::GnomadGenome3_1) {

      auto snp_results = non_const_population->inSituFilter(SNPFilter());
      ExecEnv::log().info("Population: {}, total variants: {}, 'SNP' variants: {}",
                          non_const_population->populationId(), snp_results.first, snp_results.second);
      auto unique_results = non_const_population->inSituFilter(UniqueUnphasedFilter());
      ExecEnv::log().info("Population: {}, total variants: {}, 'Unique' variants: {}",
                          non_const_population->populationId(), unique_results.first, unique_results.second);

      ExecEnv::log().info("Check Structure of Population: {}", non_const_population->populationId());

      for (auto const& [genome_id, genome_ptr] : non_const_population->getMap()) {

        for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

          for (auto const& [offset, offset_ptr] : contig_ptr->getMap()) {

            auto const& offsetVariants = offset_ptr->getVariantArray();

            size_t offset_size = offsetVariants.size();
            if (offset_size > 2) {

              for (auto const& variant_ptr : offsetVariants) {

                ExecEnv::log().info("Unexpected Count: {}, Genome: {}, {}", offset_size, genome_id, variant_ptr->output(',', VariantOutputIndex::START_0_BASED, true));

              }

            } else if (offset_size == 2){

              if (offsetVariants.front()->phaseId() == offsetVariants.back()->phaseId() and offsetVariants.front()->phaseId() != VariantPhase::UNPHASED) {

                for (auto const& variant_ptr : offsetVariants) {

                  ExecEnv::log().info("Incorrect Phasing Genome: {}, {}", genome_id, variant_ptr->output(',', VariantOutputIndex::START_0_BASED, true));

                }

              }

            }

          }

        }

      }

    } else {

      auto unique_variants = non_const_population->filterVariants(UniqueUnphasedFilter());
      ExecEnv::log().info("Population: {}, 'Unique' variants: {}",
                          non_const_population->populationId(), unique_variants->variantCount());


    }

    auto pass_results = non_const_population->inSituFilter(PassFilter());
    ExecEnv::log().info("Population: {}, total variants: {}, 'Pass' variants: {}",
                        non_const_population->populationId(), pass_results.first, pass_results.second);

    ExecEnv::log().info("Population: {}, check total variants: {}",
                        non_const_population->populationId(), non_const_population->variantCount());

    auto snp_filter_count = non_const_population->inSituFilter(SNPFilter());
    ExecEnv::log().info("Population: {}, total variants: {}, SNP variants: {}",
                        non_const_population->populationId(), snp_filter_count.first, snp_filter_count.second);

    ExecEnv::log().info("Population: {}, check total variants: {}",
                        non_const_population->populationId(), non_const_population->variantCount());

    // Delete the population by applying the FalseFilter to the population, this appears to be faster than unwinding smart pointers.
    auto false_filter_count = non_const_population->inSituFilter(FalseFilter());
    ExecEnv::log().info("Population: {}, total variants: {}, False variants: {}",
                        non_const_population->populationId(), false_filter_count.first, false_filter_count.second);

    ExecEnv::log().info("Population: {}, check total variants: {}",
                        non_const_population->populationId(), non_const_population->variantCount());

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::VerifyAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::VerifyAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
