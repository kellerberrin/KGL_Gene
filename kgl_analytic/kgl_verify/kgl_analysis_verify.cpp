//
// Created by kellerberrin on 7/11/20.
//

#include "kgl_analysis_verify.h"
#include "kgl_variant_filter.h"

#include <chrono>
#include <thread>

namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::VerifyAnalysis::initializeAnalysis( const std::string& work_directory,
                                              const ActiveParameterList& named_parameters,
                                              const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident);

  }

  for (auto const& genome_resource_ptr : resource_ptr->getResources(ResourceProperties::GENOME_RESOURCE_ID_)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());

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

      auto pass_results = non_const_population->inSituFilter(PassFilter());
      ExecEnv::log().info("Population: {}, pre filter count: {}, 'Pass' variants: {}",
                          non_const_population->populationId(), pass_results.first, pass_results.second);

      auto diploid_results = non_const_population->inSituFilter(DiploidFilter());
      ExecEnv::log().info("Population: {}, pre filter count: {}, 'Diploid' variants: {}",
                          non_const_population->populationId(), diploid_results.first, diploid_results.second);

      ExecEnv::log().info("Check Structure of Population: {}", non_const_population->populationId());

      bool check_diploid{true};
      if (check_diploid) {

        size_t incorrect_phasing{0};
        size_t unexpected_count{0};
        bool print_variants{true};
        std::map<ContigOffset_t, size_t> allele_offset_map;
        for (auto const& [genome_id, genome_ptr] : non_const_population->getMap()) {

          for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

            for (auto const& [offset, offset_ptr] : contig_ptr->getMap()) {

              auto const& offsetVariants = offset_ptr->getVariantArray();

              for (auto const& variant : offsetVariants) {

                auto result = allele_offset_map.find(variant->alleleOffset());
                if (result != allele_offset_map.end()) {

                  auto& [allele_offset, count] = *result;
                  ++count;

                } else {

                  auto insert_result = allele_offset_map.try_emplace(variant->alleleOffset(), 1);
                  const auto& [iter, bool_result] = insert_result;
                  if (not bool_result) {

                    ExecEnv::log().error("Could not insert into allele count map");

                  }

                }

              }

              size_t offset_size = offsetVariants.size();
              if (offset_size > 2) {

                unexpected_count += offset_size;
                if (print_variants) {

                  for (auto const& variant_ptr : offsetVariants) {

                    ExecEnv::log().info("Unexpected Count: {}, Genome: {}, {}, Hash: {}",
                                        offset_size,
                                        genome_id,
                                        variant_ptr->output(',', VariantOutputIndex::START_0_BASED, true),
                                        variant_ptr->HGVS());

                  }

                }

              } else if (offset_size == 2){

                if (offsetVariants.front()->phaseId() == offsetVariants.back()->phaseId() and offsetVariants.front()->phaseId() != VariantPhase::UNPHASED) {

                  incorrect_phasing += offset_size;
                  if (print_variants) {

                    for (auto const& variant_ptr : offsetVariants) {

                      ExecEnv::log().info("Incorrect Phasing Genome: {}, {}, Hash: {}",
                                          genome_id,
                                          variant_ptr->output(',', VariantOutputIndex::START_0_BASED, true),
                                          variant_ptr->HGVS());

                    }

                  }

                }

              }

            }

          }

        }

        ExecEnv::log().info("Population: {} verify check found incorrect phased variants: {}, unexpected variant count (should be 2): {}",
                            non_const_population->populationId(), incorrect_phasing, unexpected_count);

        for (auto const& [allele_offset, count] : allele_offset_map) {

          ExecEnv::log().info("Population: {} allele offset: {}, count: {}", non_const_population->populationId(), allele_offset, count);

        }

      }

    } else {

//      auto unique_variants_ = non_const_population->filterVariants(UniqueUnphasedFilter());
//      ExecEnv::log().info("Population: {}, 'Unique' variants: {}",
//                          non_const_population->populationId(), unique_variants_->variantCount());


    }


    std::pair<size_t, size_t> mem_pair = Utility::process_mem_usage2(); // pair.first is process vm_usage, pair.second is resident memory set.
    ExecEnv::log().info("Population: {}, check total variants: {}, Variant Objects:{}, Data Blocks:{}, VM Usage: {}, Resident Memory: {}",
                        non_const_population->populationId(), non_const_population->variantCount(), Variant::objectCount(),
                        DataMemoryBlock::objectCount(), mem_pair.first, mem_pair.second);

    auto empty_count = non_const_population->inSituFilter(FalseFilter());

    AuditMemory::trimFreeStore();

    const size_t seconds = 60;
    ExecEnv::log().info("*********** Sleep for seconds: {} ******************", seconds);
    std::chrono::seconds timespan(seconds); // or whatever
    std::this_thread::sleep_for(timespan);

    mem_pair = Utility::process_mem_usage2(); // pair.first is process vm_usage, pair.second is resident memory set.
    ExecEnv::log().info("Population: {}, Before variants: {}, After variants: {}, Variant Objects: {}, Data Blocks:{}, VM Usage: {}, Resident Memory: {}",
                        non_const_population->populationId(), empty_count.first, empty_count.second, Variant::objectCount(),
                        DataMemoryBlock::objectCount(), mem_pair.first, mem_pair.second);

    ExecEnv::log().info("Population: {}, Allocated Bytes: {}, Deallocated Bytes: {}, Allocations: {}, Deallocations: {}, Not max byte aligned: {}, Max aligned additional bytes: {}",
                        non_const_population->populationId(), AuditMemory::allocatedBytes(), AuditMemory::deallocatedBytes(),
                        AuditMemory::allocations(), AuditMemory::deallocations(), (AuditMemory::allocations() - AuditMemory::countMaxAlign()),
                        AuditMemory::additionalAlignBytes());


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

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
