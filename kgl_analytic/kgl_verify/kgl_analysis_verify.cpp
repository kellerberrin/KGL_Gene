//
// Created by kellerberrin on 7/11/20.
//

#include "kgl_analysis_verify.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::VerifyAnalysis::initializeAnalysis( const std::string& work_directory,
                                              const RuntimeParameterMap& named_parameters,
                                              std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

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

    // If we have a mono genome file or phased diploid then we do not expect to see duplicate variants.
    // Note that a phased diploid population will have variants distinguished by phasing.
    if (file_characteristic.data_structure == DataStructureEnum::UnphasedMonoGenome
        or file_characteristic.data_structure == DataStructureEnum::DiploidPhased) {

      ExecEnv::log().info("Verifying Population: {} for duplicate variants", population->populationId());
      VerifyDuplicates verify_duplicates;
      population->processAll(verify_duplicates, &VerifyDuplicates::verifyVariant);
      ExecEnv::log().info("Completed Verifying population for duplicate variants, duplicates: {}",
                          verify_duplicates.duplicateCount());

    }

    // Apply some filters inSitu, note this modifies the population structure so we need to cast away const-ness
    auto non_const_population = std::const_pointer_cast<PopulationDB>(population);

    auto pass_results = non_const_population->inSituFilter(PassFilter());
    ExecEnv::log().info("Population: {}, total variants: {}, 'Pass' variants: {}",
                        non_const_population->populationId(), pass_results.first, pass_results.second);

    auto snp_filter_count = non_const_population->inSituFilter(SNPFilter());
    ExecEnv::log().info("Population: {}, total variants: {}, SNP variants: {}",
                        non_const_population->populationId(), snp_filter_count.first, snp_filter_count.second);

    // Delete the population by applying the FalseFilter to the population, this appears to be faster than unwinding smart pointers.
    auto false_filter_count = non_const_population->inSituFilter(FalseFilter());
    ExecEnv::log().info("Population: {}, total variants: {}, False variants: {}",
                        non_const_population->populationId(), false_filter_count.first, false_filter_count.second);

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

bool kgl::VerifyDuplicates::verifyVariant(const std::shared_ptr<const Variant> variant_ptr) {

  if (not previous_variant_ptr_) {

    previous_variant_ptr_ = variant_ptr;
    offset_variants_.push_back(variant_ptr);
    return true;

  }

  if (previous_variant_ptr_->offset() != variant_ptr->offset()
      or previous_variant_ptr_->contigId() != variant_ptr->contigId()) {

    verifyDuplicates();
    offset_variants_.clear();

  }

  previous_variant_ptr_ = variant_ptr;
  offset_variants_.push_back(variant_ptr);

  return true;

}

bool kgl::VerifyDuplicates::verifyDuplicates() {

  const size_t offset_size = offset_variants_.size();

  for (size_t idx1 = 0; idx1 < offset_size; ++idx1) {

    for (size_t idx2 = idx1 + 1; idx2 < offset_size; ++idx2) {

      if (offset_variants_[idx1]->analogous(*offset_variants_[idx2])) {

        ExecEnv::log().info("VerifyDuplicates::verifyDuplicates; duplicate variants *************");
        ExecEnv::log().info("VerifyDuplicates::verifyDuplicates; duplicate variant : {}",
                            offset_variants_[idx1]->output(',', VariantOutputIndex::START_0_BASED, true));
        ExecEnv::log().info("VerifyDuplicates::verifyDuplicates; duplicate variant : {}",
                            offset_variants_[idx2]->output(',', VariantOutputIndex::START_0_BASED, true));

        ++duplicate_count_;

      }

    }

  }

  return true;

}
