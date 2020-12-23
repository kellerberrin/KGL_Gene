//
// Created by kellerberrin on 7/11/20.
//

#include "kgl_analysis_verify.h"
#include "kgl_variant_db_phased.h"
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
bool kgl::VerifyAnalysis::fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_ptr) {

  ExecEnv::log().info("VCF File Read for Analysis Id: {} called with Variant Population", ident());

  // List all the Info fields to remind us what's available.
  std::shared_ptr<const InfoEvidenceHeader> evidence_header_ptr;

  if (data_ptr->dataType() == DataTypeEnum::DiploidPopulation) {

    std::shared_ptr<const DiploidPopulation> diploid_population = std::dynamic_pointer_cast<const DiploidPopulation>(data_ptr);
    if (diploid_population) {

      // Only check variants that have passed all VCF filters.
      std::shared_ptr<const DiploidPopulation> diploid_passed = diploid_population->filterVariants(PassFilter());
      size_t count_passed = diploid_passed->variantCount();
      size_t count_all = diploid_population->variantCount();
      ExecEnv::log().info("Diploid Population: {}, total variants: {}, variants 'Pass' VCF filters: {}",
                          diploid_passed->populationId(), count_all, count_passed);

    }

  }

  if (data_ptr->dataType() == DataTypeEnum::UnphasedPopulation) {

    // List all the Info fields to remind us what's available.
    // Superclass the population
    std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(data_ptr);
    if (unphased_population) {

      // Only check variants that have passed all VCF filters.
      std::shared_ptr<const UnphasedPopulation> unphased_passed = unphased_population->mtFilterVariants(PassFilter());
      size_t count_passed = unphased_passed->variantCount();
      size_t count_all = unphased_population->variantCount();
      ExecEnv::log().info("Unphased Population: {}, total variants: {}, variants 'Pass' VCF filters: {}",
                          unphased_passed->populationId(), count_all, count_passed);

//      ExecEnv::log().info("Verifying Unphased 'Pass' Population: {} for duplicate variants", unphased_passed->populationId());
//      VerifyDuplicates verify_duplicates;
//      unphased_passed->processAll(verify_duplicates, &VerifyDuplicates::verifyVariant);
//      ExecEnv::log().info("Completed Verifying Unphased 'Pass' population for duplicate variants, duplicates: {}",
//                          verify_duplicates.duplicateCount());

    }

  }

  if (data_ptr->dataType() == DataTypeEnum::HaploidPopulation) {

    // List all the Info fields to remind us what's available.
    // Superclass the population
    std::shared_ptr<const UnphasedPopulation> haploid_population = std::dynamic_pointer_cast<const UnphasedPopulation>(data_ptr);
    if (haploid_population) {

      // Only check variants that have passed all VCF filters.
      std::shared_ptr<const UnphasedPopulation> haploid_passed = haploid_population->filterVariants(PassFilter());
      size_t count_passed = haploid_passed->variantCount();
      size_t count_all = haploid_population->variantCount();
      ExecEnv::log().info("Haploid Population: {}, total variants: {}, variants 'Pass' VCF filters: {}",
                          haploid_passed->populationId(), count_all, count_passed);

      ExecEnv::log().info("Verifying Haploid 'Pass' Population: {} for duplicate variants", haploid_passed->populationId());
      VerifyDuplicates verify_duplicates;
      haploid_passed->processAll(verify_duplicates, &VerifyDuplicates::verifyVariant);
      ExecEnv::log().info("Completed Verifying Haploid 'Pass' population for duplicate variants, duplicates: {}",
                          verify_duplicates.duplicateCount());

    }


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
