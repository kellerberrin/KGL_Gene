//
// Created by kellerberrin on 4/5/20.
//

#include "kgl_analysis_null.h"
#include "kgl_variant_db_phased.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::NullAnalysis::initializeAnalysis(const std::string& work_directory,
                                              const RuntimeParameterMap& named_parameters,
                                              std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::fileReadAnalysis(std::shared_ptr<const DataObjectBase> population_base) {

  ExecEnv::log().info("Default VCF File Read for Analysis Id: {} called with Variant Population", ident());

  // Superclass the population
  std::shared_ptr<const DiploidPopulation> diploid_population = std::dynamic_pointer_cast<const DiploidPopulation>(population_base);
  std::shared_ptr<const UnphasedPopulation> unphased_population = std::dynamic_pointer_cast<const UnphasedPopulation>(population_base);

  // List all the Info fields to remind us what's available.
  std::shared_ptr<const InfoEvidenceHeader> evidence_header_ptr;

  if (diploid_population) {

   std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = diploid_population->getVCFInfoEvidenceHeader();

    if (not info_header_opt) {

      ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis could not get Info Field Header from VCF diploid population: {}",
                          diploid_population->populationId());
      return false;

    }

    evidence_header_ptr = info_header_opt.value();

  }

  if (unphased_population) {

    // List all the Info fields to remind us what's available.
    std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = unphased_population->getVCFInfoEvidenceHeader();

    if (not info_header_opt) {

      ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis could not get Info Field Header from VCF unphased population: {}",
                          unphased_population->populationId());
      return false;

    }

    evidence_header_ptr = info_header_opt.value();

  }

  if (not evidence_header_ptr) {

    ExecEnv::log().error("Analysis: {},  expected a Diploid or Unphased Population", ident());
    return false;

  }

  ExecEnv::log().info("Analysis Id: {}, VCF File vcf_population: {}, Obtained header for: {} Info Fields (listed below)",
                      ident(), evidence_header_ptr->getConstMap().size(), population_base->Id());

  // List the available INFO fields
  for (auto const&[ident, field_item] : evidence_header_ptr->getConstMap()) {

    ExecEnv::log().info("Field Id: {}, Type: {}, Number: {}, Description: {}",
                        ident, field_item.infoVCF().type, field_item.infoVCF().number,
                        field_item.infoVCF().description);

  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::NullAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::NullAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}

