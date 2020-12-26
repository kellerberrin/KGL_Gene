//
// Created by kellerberrin on 4/6/20.
//

#include "kgl_analysis_info_filter.h"

#include <fstream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///


namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kgl::InfoFilterAnalysis::initializeAnalysis( const std::string& work_directory,
                                                  const RuntimeParameterMap& named_parameters,
                                                  std::shared_ptr<const GenomeCollection> reference_genomes) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

  }

  if (not getParameters(work_directory, named_parameters)) {

    return false;

  }

// Clear the data file.
  std::ofstream outfile;
  outfile.open(output_file_name_, std::ofstream::out | std::ofstream::trunc);
  outfile.close();

  return true;

}

// Perform the genetic analysis per iteration.
bool kgl::InfoFilterAnalysis::fileReadAnalysis(std::shared_ptr<const DataObjectBase> data_ptr) {

  // Analysis currently only defined for Gnomad 2.1
  auto file_characteristic = data_ptr->dataCharacteristic();
  if (file_characteristic.data_source != DataSourceEnum::Gnomad2_1) {

    ExecEnv::log().warn("Analysis: {}, expected a Gnomad 2.1 Population for file: {}", ident(), data_ptr->fileId());
    return true;

  }

  // Superclass the population
  std::shared_ptr<const PopulationVariant> vcf_population = std::dynamic_pointer_cast<const PopulationVariant>(data_ptr);

  if (not vcf_population) {

    ExecEnv::log().info("Analysis: {}, expected a Population for file: {}", ident(), data_ptr->fileId());
    return false;

  }

  // Pre-filter variants for quality, using the VQSLOD and rf_tp_probability fields.
  filtered_vcf_population_ = qualityFilter(vcf_population);

  std::ofstream outfile(output_file_name_,  std::ofstream::out | std::ofstream::app);

  if (not outfile.good()) {

    ExecEnv::log().error("InfoFilterAnalysis::finalizeAnalysis; could not open results file: {}", output_file_name_);
    return false;

  }


  // Perform the chromosome analysis.
  bool result = performAnalysis(vcf_population);

  return result;

}



// Perform the genetic analysis per iteration.
bool kgl::InfoFilterAnalysis::iterationAnalysis() {

// No operations for an iteration.

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kgl::InfoFilterAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}



bool kgl::InfoFilterAnalysis::getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters) {

  // Get the output filename
  auto result = named_parameters.find(OUTPUT_FILE_);

  if (result == named_parameters.end()) {

    ExecEnv::log().error("Analytic: {}; Expected Parameter: {} to be defined. {} is deactivated. Available named Parameters:", ident(), OUTPUT_FILE_, ident());
    for (auto const& [parameter_ident, parameter_value] : named_parameters) {

      ExecEnv::log().info("Analysis: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

    }

    return false;

  }

  output_file_name_ = Utility::filePath(result->second, work_directory);

  ExecEnv::log().info("Analysis: {}, initialized with output file: {}", ident(), output_file_name_);

  return true;

}

std::shared_ptr<kgl::PopulationVariant>
kgl::InfoFilterAnalysis::qualityFilter( std::shared_ptr<const PopulationVariant> vcf_population) {


  size_t unfiltered = vcf_population->variantCount();
  ExecEnv::log().info("Population: {} size: {} before filtering", vcf_population->populationId(), unfiltered);

  const std::string VQSLOD_FIELD{"VQSLOD"};
  // const double VQSLOD_LEVEL{1.775};
  const double VQSLOD_LEVEL{1.2168};
  auto vqslod_filter = InfoGEQFloatFilter(VQSLOD_FIELD, VQSLOD_LEVEL);

  std::shared_ptr<PopulationVariant> filtered_population = vcf_population->filterVariants(vqslod_filter);

  size_t filtered_VQSLOD = filtered_population->variantCount();
  double percent_filtered =  (static_cast<double>(filtered_VQSLOD) / static_cast<double>(unfiltered)) * 100.0;

  ExecEnv::log().info("Population: {} size: {} ({}%) after VQSLOD filter level: {}",
                      filtered_population->populationId(), filtered_VQSLOD, percent_filtered, VQSLOD_LEVEL);

  const std::string RANDOM_FOREST_FIELD{"rf_tp_probability"};
  const double RANDOM_FOREST_LEVEL{0.90};

  auto random_forest_filter = InfoGEQFloatFilter(RANDOM_FOREST_FIELD, RANDOM_FOREST_LEVEL);

  filtered_population = filtered_population->filterVariants(random_forest_filter);

  size_t filtered_random_forest = filtered_population->variantCount();
  percent_filtered =  (static_cast<double>(filtered_random_forest) / static_cast<double>(unfiltered)) * 100.0;

  ExecEnv::log().info("Population: {} size: {} ({}%) after VQSLOD and Random Forest filter level: {}",
                      filtered_population->populationId(), filtered_population->variantCount(), percent_filtered, RANDOM_FOREST_LEVEL);

  return filtered_population;

}


bool kgl::InfoFilterAnalysis::performAnalysis( std::shared_ptr<const kgl::PopulationVariant> vcf_population) {


  std::ofstream outfile(output_file_name_,  std::ofstream::out | std::ofstream::app);

  if (not outfile.good()) {

    ExecEnv::log().error("InfoFilterAnalysis::fileReadAnalysis; could not open results file: {}", output_file_name_);
    return false;

  }


  ExecEnv::log().info("Population: {} size after filtering: {}", vcf_population->populationId(), vcf_population->variantCount());

  //perform the analysis.
  InfoAgeAnalysis age_sub_total(" Total for Population: " + vcf_population->populationId());
  vcf_population->processAll(age_sub_total, &InfoAgeAnalysis::processVariant);
  outfile << age_sub_total;
  // Store for final total.
  age_analysis_vector_.push_back(age_sub_total);


  // Filter on the VEP Impact field.
  const std::string vep_sub_field("IMPACT");
  const std::string vep_high_impact("HIGH");
  const std::string vep_moderate_impact("MODERATE");

  auto high_impact_filter =  VepSubStringFilter(vep_sub_field, vep_high_impact);
  auto moderate_impact_filter = VepSubStringFilter(vep_sub_field, vep_moderate_impact);
  auto vep_impact_filter = OrFilter(high_impact_filter, moderate_impact_filter);
  // Analyze high and moderate impact variants.
  analyzeFilteredPopulation(vep_impact_filter, vcf_population, outfile);

  // Filter by age deciles.
  filterByAge(vcf_population, outfile);

  // Filter on AF fields.
  const std::vector<double> AF_values{0.001, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.90};
  const std::vector<double> inbreeding_values{ -0.9, -0.5, -0.2, -0.1, 0.0, 0.1, 0.20, 0.50, 0.90};

  analyzeField("InbreedingCoeff", inbreeding_values, vcf_population, outfile);
  analyzeField("AF", AF_values, vcf_population, outfile);
  analyzeField("AF_male", AF_values, vcf_population, outfile);
  analyzeField("AF_female", AF_values, vcf_population, outfile);
  analyzeField("AF_afr", AF_values, vcf_population, outfile);
  analyzeField("AF_amr", AF_values, vcf_population, outfile);
  analyzeField("AF_asj", AF_values, vcf_population, outfile);
  analyzeField("AF_eas", AF_values, vcf_population, outfile);
  analyzeField("AC_eas_jpn", AF_values, vcf_population, outfile);
  analyzeField("AC_eas_kor", AF_values, vcf_population, outfile);
  analyzeField("AC_eas_oea", AF_values, vcf_population, outfile);
  analyzeField("AF_fin", AF_values, vcf_population, outfile);
  analyzeField("AF_nfe", AF_values, vcf_population, outfile);
  analyzeField("AF_nfe_est", AF_values, vcf_population, outfile);
  analyzeField("AC_nfe_bgr", AF_values, vcf_population, outfile);
  analyzeField("AF_nfe_nwe", AF_values, vcf_population, outfile);
  analyzeField("AF_nfe_onf", AF_values, vcf_population, outfile);
  analyzeField("AF_nfe_seu", AF_values, vcf_population, outfile);
  analyzeField("AC_nfe_swe", AF_values, vcf_population, outfile);
  analyzeField("AF_oth", AF_values, vcf_population, outfile);

  return true;

}


void kgl::InfoFilterAnalysis::listAvailableInfoFields(std::shared_ptr<const PopulationVariant> vcf_population) {

  // Investigate vep field values.
  InfoEvidenceAnalysis::vepSubFieldValues("Consequence", vcf_population);
  InfoEvidenceAnalysis::vepSubFieldValues("IMPACT", vcf_population);
  InfoEvidenceAnalysis::vepSubFieldValues("Feature_type", vcf_population);
  InfoEvidenceAnalysis::vepSubFieldValues("BIOTYPE", vcf_population);

  // List all the Info fields to remind us what's available.
  std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = vcf_population->getVCFInfoEvidenceHeader();

  if (not info_header_opt) {

    ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis could not get Info Field Header from VCF vcf_population: {}",
                        vcf_population->populationId());
    return;

  }

  ExecEnv::log().info("Analysis Id: {}, VCF File vcf_population: {}, Obtained header for: {} Info Fields (listed below)",
                      ident(), info_header_opt.value()->getConstMap().size(), vcf_population->populationId());

  // List the available fields for ease of use.
  for (auto const&[ident, field_item] : info_header_opt.value()->getConstMap()) {

    ExecEnv::log().info("Field Id: {}, Type: {}, Number: {}, Description: {}",
                        ident, field_item.infoVCF().type, field_item.infoVCF().number,
                        field_item.infoVCF().description);

  }


}


void kgl::InfoFilterAnalysis::filterByAge( std::shared_ptr<const PopulationVariant> vcf_population,
                                           std::ostream& result_file) {

  ExecEnv::log().info("Perform Homo/Hetero analysis with average variant ages filtered to 10 deciles.");

  AgeSortedMap age_sorted_map;

  vcf_population->processAll(age_sorted_map, &AgeSortedMap::processVariant);
  // Save for the final analysis.
  vcf_population->processAll(age_sorted_map_, &AgeSortedMap::processVariant);

  result_file << age_sorted_map;

}


void kgl::InfoFilterAnalysis::analyzeField( const std::string& info_field_ident,
                                            const std::vector<double>& field_values,
                                            std::shared_ptr<const PopulationVariant> vcf_population,
                                            std::ostream& result_file) {

  analyzeFilteredPopulation(NotFilter(InfoGEQFloatFilter(info_field_ident, field_values.front())), vcf_population, result_file);

  for (auto value : field_values) {

    analyzeFilteredPopulation(InfoGEQFloatFilter(info_field_ident, value), vcf_population, result_file);

  }

}



void kgl::InfoFilterAnalysis::analyzeFilteredPopulation( const VariantFilter& filter,
                                                         std::shared_ptr<const PopulationVariant> vcf_population,
                                                         std::ostream& result_file) {

  // Tag with filter and population
  std::string title = "Filter: " + filter.filterName() + ", Population: " + vcf_population->populationId();
  ExecEnv::log().info("Analysis Package: {}, executing age analysis: {}", ident(), title);
  InfoAgeAnalysis age_analysis(title);
  // Filter the variant population
  std::shared_ptr<const PopulationVariant> filtered_population = vcf_population->filterVariants(filter);
  // Gather the age profile.
  filtered_population->processAll(age_analysis, &InfoAgeAnalysis::processVariant);
  // Write results
  result_file << age_analysis;

}

