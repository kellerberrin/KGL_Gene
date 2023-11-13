//
// Created by kellerberrin on 4/6/20.
//

#include "kga_analysis_info_filter.h"

#include <fstream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///


namespace kga = kellerberrin::genome::analysis;
namespace kgl = kellerberrin::genome;


// Setup the analytics to process VCF data.
bool kga::InfoFilterAnalysis::initializeAnalysis( const std::string& work_directory,
                                                  const ActiveParameterList& named_parameters,
                                                  const std::shared_ptr<const AnalysisResources>& resource_ptr) {

  ExecEnv::log().info("Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_map] : named_parameters.getMap()) {

    ExecEnv::log().info("Initialize Analysis Id: {}, initialized with parameter: {}", ident(), parameter_ident);

  }

  for (auto const& genome_resource_ptr : resource_ptr->getResources(ResourceProperties::GENOME_RESOURCE_ID_)) {

    auto genome_ptr = std::dynamic_pointer_cast<const GenomeReference>(genome_resource_ptr);
    ExecEnv::log().info("Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome_ptr->genomeId());

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
bool kga::InfoFilterAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> data_ptr) {

  // Analysis currently only defined for Gnomad 2.1
  auto file_characteristic = data_ptr->dataCharacteristic();
  if ( file_characteristic.data_source != DataSourceEnum::Gnomad2_1
      and file_characteristic.data_source != DataSourceEnum::Gnomad3_1
      and file_characteristic.data_source != DataSourceEnum::GnomadExomes2_1
      and file_characteristic.data_source != DataSourceEnum::GnomadExomes3_1) {

    ExecEnv::log().warn("Analysis: {}, expected a Gnomad 2.1 Population for file: {}", ident(), data_ptr->fileId());
    return true;

  }

  // Superclass the population
  std::shared_ptr<const PopulationDB> vcf_population = std::dynamic_pointer_cast<const PopulationDB>(data_ptr);

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
bool kga::InfoFilterAnalysis::iterationAnalysis() {

// No operations for an iteration.

  return true;

}

// All VCF data has been presented, finalize analysis and write results.
bool kga::InfoFilterAnalysis::finalizeAnalysis() {

  ExecEnv::log().info("Finalize Analysis called for Analysis Id: {}", ident());

  return true;

}


bool kga::InfoFilterAnalysis::getParameters( const std::string& work_directory,
                                             const ActiveParameterList& named_parameters) {


  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    if (block_vector.size() != 1) {

      ExecEnv::log().error("InfoFilterAnalysis::getParameters; parameter block: {} vector size: {}, expected size = 1",
                           block_name, block_vector.size());
      return false;

    }

    ExecEnv::log().info("Analysis: {} parsing parameter block: {}", ident(), block_name);

    for (auto const& xml_vector : block_vector) {

      auto output_opt = xml_vector.getString(OUTPUT_FILE_);
      if (output_opt) {

        output_file_name_ = output_opt.value().front() + std::string(OUTPUT_FILE_EXT_);
        output_file_name_ = Utility::filePath(output_file_name_, work_directory);


      } else {

        ExecEnv::log().error("InfoFilterAnalysis::getParameters; bad value for parameter: {}", OUTPUT_FILE_);
        return false;

      }

    }

  }

  return true;

}

std::shared_ptr<kgl::PopulationDB>
kga::InfoFilterAnalysis::qualityFilter( std::shared_ptr<const PopulationDB> vcf_population) {


  size_t unfiltered = vcf_population->variantCount();
  ExecEnv::log().info("Population: {} size: {} before filtering", vcf_population->populationId(), unfiltered);

  const std::string VQSLOD_FIELD{"VQSLOD"};
  // const double VQSLOD_LEVEL{1.775};
  const double VQSLOD_LEVEL{1.2168};
  auto vqslod_filter = InfoFilter<double, false>(VQSLOD_FIELD, [VQSLOD_LEVEL](double compare) ->bool { return compare >= VQSLOD_LEVEL; });

  std::shared_ptr<PopulationDB> filtered_population = vcf_population->viewFilter(vqslod_filter);

  size_t filtered_VQSLOD = filtered_population->variantCount();
  double percent_filtered =  (static_cast<double>(filtered_VQSLOD) / static_cast<double>(unfiltered)) * 100.0;

  ExecEnv::log().info("Population: {} size: {} ({}%) after VQSLOD filter level: {}",
                      filtered_population->populationId(), filtered_VQSLOD, percent_filtered, VQSLOD_LEVEL);

  const std::string RANDOM_FOREST_FIELD{"rf_tp_probability"};
  const double RANDOM_FOREST_LEVEL{0.90};

  auto random_forest_filter = InfoFilter<double, false>(RANDOM_FOREST_FIELD, [RANDOM_FOREST_LEVEL](double compare)->bool { return compare >= RANDOM_FOREST_LEVEL; });

  filtered_population = filtered_population->viewFilter(random_forest_filter);

  size_t filtered_random_forest = filtered_population->variantCount();
  percent_filtered =  (static_cast<double>(filtered_random_forest) / static_cast<double>(unfiltered)) * 100.0;

  ExecEnv::log().info("Population: {} size: {} ({}%) after VQSLOD and Random Forest filter level: {}",
                      filtered_population->populationId(), filtered_population->variantCount(), percent_filtered, RANDOM_FOREST_LEVEL);

  return filtered_population;

}


bool kga::InfoFilterAnalysis::performAnalysis( std::shared_ptr<const kgl::PopulationDB> vcf_population) {


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

  auto high_impact_filter =  VepSubStringFilter<false>(vep_sub_field, vep_high_impact);
  auto moderate_impact_filter = VepSubStringFilter<false>(vep_sub_field, vep_moderate_impact);
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


void kga::InfoFilterAnalysis::listAvailableInfoFields(std::shared_ptr<const PopulationDB> vcf_population) {

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


void kga::InfoFilterAnalysis::filterByAge( std::shared_ptr<const PopulationDB> vcf_population,
                                           std::ostream& result_file) {

  ExecEnv::log().info("Perform Homo/Hetero analysis with average variant ages filtered to 10 deciles.");

  AgeSortedMap age_sorted_map;

  vcf_population->processAll(age_sorted_map, &AgeSortedMap::processVariant);
  // Save for the final analysis.
  vcf_population->processAll(age_sorted_map_, &AgeSortedMap::processVariant);

  result_file << age_sorted_map;

}


void kga::InfoFilterAnalysis::analyzeField( const std::string& info_field_ident,
                                            const std::vector<double>& field_values,
                                            std::shared_ptr<const PopulationDB> vcf_population,
                                            std::ostream& result_file) {

  if (field_values.empty()) {

    ExecEnv::log().error("nfoFilterAnalysis::analyzeField; Unaexpected empty value vector");
    return;

  }

  // Comparison lambda returns true if compare is less than the front value.
  auto float_front_lambda = [value=field_values.front()](double compare)->bool { return compare < value; };
  // Comparison filter returns false if the info_field_ident does not exist.
  InfoFilter<double, false> compare_front_filter(info_field_ident, float_front_lambda);

  analyzeFilteredPopulation(compare_front_filter, vcf_population, result_file);

  for (auto value : field_values) {

    // Comparison lambda returns true if compare is greater or equal than the iterated field value.
    auto float_value_lambda = [value](double compare)->bool { return compare >= value; };

    // Comparison filter returns false if the info_field_ident does not exist.
    InfoFilter<double, false> compare_value_filter(info_field_ident, float_value_lambda);

    analyzeFilteredPopulation(compare_value_filter, vcf_population, result_file);

  }

}



void kga::InfoFilterAnalysis::analyzeFilteredPopulation( const BaseFilter& filter,
                                                         std::shared_ptr<const PopulationDB> vcf_population,
                                                         std::ostream& result_file) {

  // Tag with filter and population
  std::string title = "Filter: " + filter.filterName() + ", Population: " + vcf_population->populationId();
  ExecEnv::log().info("Analysis Package: {}, executing age analysis: {}", ident(), title);
  InfoAgeAnalysis age_analysis(title);
  // Filter the variant population
  auto filtered_population_ptr = vcf_population->viewFilter(filter);
  // Gather the age profile.
  filtered_population_ptr->processAll(age_analysis, &InfoAgeAnalysis::processVariant);
  // Write results
  result_file << age_analysis;

}

