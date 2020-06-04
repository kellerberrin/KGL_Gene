//
// Created by kellerberrin on 4/6/20.
//

#include "kgl_analysis_info_filter.h"


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::InfoFilterAnalysis::initializeAnalysis( const std::string& work_directory,
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
bool kgl::InfoFilterAnalysis::fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_population) {

  ExecEnv::log().info("Default VCF File Read for Analysis Id: {} called with Variant Population: {}, Variant Count: {}",
                      ident(), vcf_population->populationId(), vcf_population->variantCount());


  std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = vcf_population->getVCFInfoEvidenceHeader();

  if (info_header_opt) {

    ExecEnv::log().info("Analysis Id: {}, VCF File vcf_population: {}, Obtained header for: {} Info Fields (listed below)",
                        ident(), info_header_opt.value()->getConstMap().size(), vcf_population->populationId());

    for (auto const& [ident, field_item] : info_header_opt.value()->getConstMap()) {

      ExecEnv::log().info( "Field Id: {}, Type: {}, Number: {}, Description: {}",
                           ident, field_item.infoVCF().type, field_item.infoVCF().number, field_item.infoVCF().description);

    }

    std::string filter_field_ident("AF");

    std::optional<const InfoSubscribedField> filter_field = info_header_opt.value()->getSubscribedField(filter_field_ident);

    if (filter_field) {


      InfoGEQFloatFilter info_filter(filter_field.value(), 0.01, false);

      ExecEnv::log().info( "Filtering VCF Population using field: {}", filter_field_ident);

      std::shared_ptr<const UnphasedPopulation> filtered_population = vcf_population->filterVariants(info_filter);

      ExecEnv::log().info( "Filter complete, Original VCF Population size: {}, Filtered Population size: {}",
                            vcf_population->variantCount(), filtered_population->variantCount());

    } else {

      ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis, filter Info Field: {} not found. Disabled.", filter_field_ident);
      return false;

    }


  } else {

    ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis could not get Info Field Header from VCF vcf_population: {}. Disabled.",
                        vcf_population->populationId());
    return false;


  }

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::InfoFilterAnalysis::iterationAnalysis() {

  ExecEnv::log().info("Default Iteration Analysis called for Analysis Id: {}", ident());

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::InfoFilterAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Default Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}

