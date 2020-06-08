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

  ExecEnv::log().info("Default Analysis Id: {} initialized with work directory: {}", ident(), work_directory);
  for (auto const& [parameter_ident, parameter_value] : named_parameters) {

    ExecEnv::log().info("Default Initialize Analysis Id: {}, initialized with parameter: {}, value: {}", ident(), parameter_ident, parameter_value);

  }

  for (auto const& genome : reference_genomes->getMap()) {

    ExecEnv::log().info("Default Initialize for Analysis Id: {} called with Reference Genome: {}", ident(), genome.first);

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
bool kgl::InfoFilterAnalysis::fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_population) {

  ExecEnv::log().info("Default VCF File Read for Analysis Id: {} called with Variant Population: {}, Variant Count: {}",
                      ident(), vcf_population->populationId(), vcf_population->variantCount());


  std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = vcf_population->getVCFInfoEvidenceHeader();

  if (info_header_opt) {

    ExecEnv::log().info("Analysis Id: {}, VCF File vcf_population: {}, Obtained header for: {} Info Fields (listed below)",
                        ident(), info_header_opt.value()->getConstMap().size(), vcf_population->populationId());

    // List the available fields ease of use.
    for (auto const& [ident, field_item] : info_header_opt.value()->getConstMap()) {

      ExecEnv::log().info( "Field Id: {}, Type: {}, Number: {}, Description: {}",
                           ident, field_item.infoVCF().type, field_item.infoVCF().number, field_item.infoVCF().description);

    }

    // Pre-filter variants for quality, uses the VQSLOD and rf_tp_probability fields.
    std::optional<std::shared_ptr<const kgl::UnphasedPopulation>> vcf_population_opt = qualityFilter(vcf_population);

    if (not vcf_population_opt) {

      ExecEnv::log().error( "InfoFilterAnalysis::fileReadAnalysis; problem with variant quality filter for population: {}",
                            vcf_population->populationId());
      return false;

    }

    vcf_population = vcf_population_opt.value();

    ExecEnv::log().info( "Population: {} size after filtering: {}", vcf_population->populationId(), vcf_population->variantCount());

    std::ofstream outfile(output_file_name_,  std::ofstream::out | std::ofstream::app);

    if (not outfile.good()) {

      ExecEnv::log().error("InfoFilterAnalysis::fileReadAnalysis; could not open results file: {}", output_file_name_);
      return false;

    }

/*
AF_afr, Type: Float, Number: A, Description: Alternate allele frequency in samples of African-American ancestry
AF_amr, Type: Float, Number: A, Description: Alternate allele frequency in samples of Latino ancestry
AF_asj, Type: Float, Number: A, Description: Alternate allele frequency in samples of Ashkenazi Jewish ancestry
AF_eas, Type: Float, Number: A, Description: Alternate allele frequency in samples of East Asian ancestry
AF_female, Type: Float, Number: A, Description: Alternate allele frequency in female samples
AF_fin, Type: Float, Number: A, Description: Alternate allele frequency in samples of Finnish ancestry
AF_male, Type: Float, Number: A, Description: Alternate allele frequency in male samples
AF_nfe, Type: Float, Number: A, Description: Alternate allele frequency in samples of non-Finnish European ancestry
AF_nfe_est, Type: Float, Number: A, Description: Alternate allele frequency in samples of Estonian ancestry
AF_nfe_nwe, Type: Float, Number: A, Description: Alternate allele frequency in samples of North-Western European ancestry
AF_nfe_onf, Type: Float, Number: A, Description: Alternate allele frequency in samples of non-Finnish but otherwise indeterminate European ancestry
AF_nfe_seu, Type: Float, Number: A, Description: Alternate allele frequency in samples of Southern European ancestry
AF_oth, Type: Float, Number: A, Description: Alternate allele frequency in samples of uncertain ancestry
*/

    InfoAgeAnalysis age_sub_total(" Total for Population: " + vcf_population->populationId());
    vcf_population->processAll(age_sub_total);
    outfile << age_sub_total;
    // Store for final total.
    age_analysis_vector_.push_back(age_sub_total);

    // Filter on AF fields.
    analyzeField("InbreedingCoeff", vcf_population, outfile);
    analyzeField("AF", vcf_population, outfile);
    analyzeField("AF_male", vcf_population, outfile);
    analyzeField("AF_female", vcf_population, outfile);
    analyzeField("AF_afr", vcf_population, outfile);
    analyzeField("AF_amr", vcf_population, outfile);
    analyzeField("AF_asj", vcf_population, outfile);
    analyzeField("AF_eas", vcf_population, outfile);
    analyzeField("AF_fin", vcf_population, outfile);
    analyzeField("AF_nfe", vcf_population, outfile);
    analyzeField("AF_nfe_est", vcf_population, outfile);
    analyzeField("AF_nfe_nwe", vcf_population, outfile);
    analyzeField("AF_nfe_onf", vcf_population, outfile);
    analyzeField("AF_nfe_seu", vcf_population, outfile);
    analyzeField("AF_oth", vcf_population, outfile);

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

  std::ofstream outfile(output_file_name_,  std::ofstream::out | std::ofstream::app);

  if (not outfile.good()) {

    ExecEnv::log().error("InfoFilterAnalysis::finalizeAnalysis; could not open results file: {}", output_file_name_);
    return false;

  }

  InfoAgeAnalysis all_age_analysis("All Populations / All Contigs");

  for (auto const& age_analysis : age_analysis_vector_) {

    all_age_analysis.addAgeAnalysis(age_analysis);

  }

  outfile << all_age_analysis << '\n';

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

std::optional<std::shared_ptr<const kgl::UnphasedPopulation>>
kgl::InfoFilterAnalysis::qualityFilter( std::shared_ptr<const UnphasedPopulation> vcf_population) {

  const std::string VQSLOD_FIELD{"VQSLOD"};
  const double VQSLOD_LEVEL{1.775};
  const std::string RANDOM_FOREST_FIELD{"rf_tp_probability"};
  const double RANDOM_FOREST_LEVEL{0.95};

  size_t unfiltered = vcf_population->variantCount();

  ExecEnv::log().info( "Population: {} size: {} before filtering", vcf_population->populationId(), unfiltered);

  std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = vcf_population->getVCFInfoEvidenceHeader();

  if (not info_header_opt) {

    return std::nullopt;

  }

  // Construct VQSLOD filter
  std::optional<const InfoSubscribedField> vqslod_field = info_header_opt.value()->getSubscribedField(VQSLOD_FIELD);

  if (not vqslod_field) {

    return std::nullopt;

  }

  InfoGEQFloatFilter vqslod_filter(vqslod_field.value(), VQSLOD_LEVEL);

  std::shared_ptr<const UnphasedPopulation> filtered_population = vcf_population->filterVariants(vqslod_filter);

  size_t filtered_VQSLOD = filtered_population->variantCount();

  double percent_filtered =  (static_cast<double>(filtered_VQSLOD) / static_cast<double>(unfiltered)) * 100.0;

  ExecEnv::log().info( "Population: {} size: {} ({}%) after VQSLOD filter level: {}",
                       filtered_population->populationId(), filtered_VQSLOD, percent_filtered, VQSLOD_LEVEL);

  // Construct Random Forest filter.
  std::optional<const InfoSubscribedField> random_forest_field = info_header_opt.value()->getSubscribedField(RANDOM_FOREST_FIELD);

  if (not random_forest_field) {

    return std::nullopt;

  }

  InfoGEQFloatFilter random_forest_filter(random_forest_field.value(), RANDOM_FOREST_LEVEL);

  filtered_population = filtered_population->filterVariants(vqslod_filter);

  size_t filtered_random_forest = filtered_population->variantCount();

  percent_filtered =  (static_cast<double>(filtered_random_forest) / static_cast<double>(unfiltered)) * 100.0;

  ExecEnv::log().info( "Population: {} size: {} ({}%) after VQSLOD and Random Forest filter level: {}",
                       filtered_population->populationId(), filtered_population->variantCount(), percent_filtered, RANDOM_FOREST_LEVEL);

  return filtered_population;

}


void kgl::InfoFilterAnalysis::analyzeField( const std::string& info_field_ident,
                                            std::shared_ptr<const UnphasedPopulation> vcf_population,
                                            std::ostream& result_file) {

  std::optional<std::shared_ptr<const InfoEvidenceHeader>> info_header_opt = vcf_population->getVCFInfoEvidenceHeader();

  if (info_header_opt) {

    std::optional<const InfoSubscribedField> filter_field = info_header_opt.value()->getSubscribedField(info_field_ident);

    if (filter_field) {

      const InfoSubscribedField& info_field = filter_field.value();

      analyzeFilteredPopulation(NotFilter(InfoGEQFloatFilter(info_field, 0.001)), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.001), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.01), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.02), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.05), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.10), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.20), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.50), vcf_population, result_file);
      analyzeFilteredPopulation(InfoGEQFloatFilter(info_field, 0.90), vcf_population, result_file);

    } else {

      ExecEnv::log().warn("InfoFilterAnalysis::fileReadAnalysis, filter Info Field: {} not found. Disabled.", info_field_ident);

    }

  }

}



void kgl::InfoFilterAnalysis::analyzeFilteredPopulation( const VariantFilter& filter,
                                                         std::shared_ptr<const UnphasedPopulation> vcf_population,
                                                         std::ostream& result_file) {

  // Tag with filter and population
  std::string title = "Filter: " + filter.filterName() + ", Population: " + vcf_population->populationId();
  ExecEnv::log().info("Analysis Package: {}, executing age analysis: {}", ident(), title);
  InfoAgeAnalysis age_analysis(title);
  // Filter the variant population
  std::shared_ptr<const UnphasedPopulation> filtered_population = vcf_population->filterVariants(filter);
  // Gather the age profile.
  filtered_population->processAll(age_analysis);
  // Write results
  result_file << age_analysis;

}

