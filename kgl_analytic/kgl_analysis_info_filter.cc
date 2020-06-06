//
// Created by kellerberrin on 4/6/20.
//

#include "kgl_analysis_info_filter.h"

#include <fstream>


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///


double kgl::InfoAgeAnalysis::processField(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> field_opt = InfoEvidenceAnalysis::getInfoData(variant_ptr, field_name);

  if (field_opt) {

    std::vector<int64_t> field_vec = InfoEvidenceAnalysis::varianttoIntegers(field_opt.value());

    if (field_vec.size() != 1) {

      ExecEnv::log().warn("InfoAgeAnalysis::processVariant, Field: {} expected vector size 1, get vector size: {}",
                          field_name, field_vec.size());

    } else {

      return static_cast<double>(field_vec.front());

    }

  } else {

    ExecEnv::log().warn("InfoAgeAnalysis::processVariant, Field: {} not available", field_name);

  }

  return 0.0;

}



std::vector<double> kgl::InfoAgeAnalysis::processBin(const std::shared_ptr<const Variant>& variant_ptr, const std::string& field_name) {


  std::optional<kgl::InfoDataVariant> bin_info_opt = InfoEvidenceAnalysis::getInfoData(variant_ptr, field_name);

  if (bin_info_opt) {

    std::vector<std::string> age_string = InfoEvidenceAnalysis::varianttoStrings(bin_info_opt.value());

    std::vector<double> age_vector = InfoEvidenceAnalysis::stringBinToFloat(age_string, AGE_BIN_SIZE_);

    return age_vector;

  } else {

    ExecEnv::log().warn("InfoAgeAnalysis::processVariant, data for age bin field: {}, not available", field_name);

  }

  return std::vector<double>(AGE_BIN_SIZE_, 0.0);

}


void kgl::InfoAgeAnalysis::processVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  ++variant_count_;

  het_under_30_ += processField(variant_ptr, HETERO_UNDER30_FIELD_);
  het_80_over_ += processField(variant_ptr, HETERO_80OVER_FIELD_);
  hom_under_30_ += processField(variant_ptr, HOMO_UNDER30_FIELD_);
  hom_80_over_ += processField(variant_ptr, HOMO_80OVER_FIELD_);

  size_t index = 0;
  for (auto& age : processBin(variant_ptr, HETERO_AGE_FIELD_)) {

    het_age_vector_[index] += age;

    ++index;

  }

  index = 0;
  for (auto& age : processBin(variant_ptr, HOMO_AGE_FIELD_)) {

    hom_age_vector_[index] += age;

    ++index;

  }

}


void kgl::InfoAgeAnalysis::addAgeAnalysis(const InfoAgeAnalysis& age_analysis) {

  size_t index = 0;
  for (auto age : age_analysis.hom_age_vector_) {

    hom_age_vector_[index] += age;
    ++index;

  }

  hom_under_30_ += age_analysis.hom_under_30_;
  hom_80_over_ += age_analysis.hom_80_over_;

  index = 0;
  for (auto age : age_analysis.het_age_vector_) {

    het_age_vector_[index] += age;
    ++index;

  }

  het_under_30_ += age_analysis.het_under_30_;
  het_80_over_ += age_analysis.het_80_over_;

  variant_count_ += age_analysis.variant_count_;

}



double kgl::InfoAgeAnalysis::sumHomozygous() const {

  double sum = hom_under_30_;
  for (auto age : hom_age_vector_) {

    sum += age;

  }

  sum += hom_80_over_;

  return sum;

}

double kgl::InfoAgeAnalysis::sumHeterozygous() const {

  double sum = het_under_30_;
  for (auto age : het_age_vector_) {

    sum += age;

  }

  sum += het_80_over_;

  return sum;

}

// Utility function writes results to a stream.
std::ostream& operator<<(std::ostream& ostream, const kellerberrin::genome::InfoAgeAnalysis& age_analysis) {

  ostream << kgl::InfoAgeAnalysis::header() << '\n';
  ostream << "hom, " << age_analysis.ageHomozygousUnder30() << ", ";
  for (auto const age : age_analysis.ageHomozygousVector()) {

    ostream << age << ", ";

  }
  ostream << age_analysis.ageHomozygous80Over() << '\n';

  double sum_hom = age_analysis.sumHomozygous() * 0.01;
  ostream << "hom%, " << (age_analysis.ageHomozygousUnder30() / sum_hom) << ", ";
  for (auto const age : age_analysis.ageHomozygousVector()) {

    ostream << (age /sum_hom) << ", ";

  }
  ostream << (age_analysis.ageHomozygous80Over() /sum_hom) << '\n';


  ostream << "het, " << age_analysis.ageHeterozygousUnder30() << ", ";
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << age << ", ";

  }
  ostream << age_analysis.ageHeterozygous80Over() << '\n';

  double sum_het = age_analysis.sumHeterozygous() * 0.01;

  ostream << "het%, " << (age_analysis.ageHeterozygousUnder30() / sum_het) << ", ";
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << (age / sum_het) << ", ";

  }
  ostream << (age_analysis.ageHeterozygous80Over() / sum_het) << '\n';

  ostream << "het/hom, " << (age_analysis.ageHeterozygousUnder30() / age_analysis.ageHomozygousUnder30()) << ", ";
  size_t index = 0;
  for (auto const age : age_analysis.ageHeterozygousVector()) {

    ostream << (age / age_analysis.ageHomozygousVector()[index]) << ", ";

    ++index;

  }
  ostream << (age_analysis.ageHeterozygous80Over() / age_analysis.ageHomozygous80Over()) << '\n';

  return ostream;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///


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

    for (auto const& [ident, field_item] : info_header_opt.value()->getConstMap()) {

      ExecEnv::log().info( "Field Id: {}, Type: {}, Number: {}, Description: {}",
                           ident, field_item.infoVCF().type, field_item.infoVCF().number, field_item.infoVCF().description);

    }

    std::string filter_field_ident("AF");

    std::optional<const InfoSubscribedField> filter_field = info_header_opt.value()->getSubscribedField(filter_field_ident);

    if (filter_field) {

      InfoAgeAnalysis age_analysis_all;
      InfoAgeAnalysis age_analysis_0_1_percent;
      InfoAgeAnalysis age_analysis_1_percent;
      InfoAgeAnalysis age_analysis_2_percent;
      InfoAgeAnalysis age_analysis_5_percent;
      InfoAgeAnalysis age_analysis_10_percent;

      vcf_population->processAll(age_analysis_all);
      InfoGEQFloatFilter info_0_1_filter(filter_field.value(), 0.001, false);
      std::shared_ptr<const UnphasedPopulation> filtered_population = vcf_population->filterVariants(info_0_1_filter);
      ExecEnv::log().info( "{}, Original VCF Population size: {}, Filtered Population size: {}",
                           info_0_1_filter.filterName(), vcf_population->variantCount(), filtered_population->variantCount());
      filtered_population->processAll(age_analysis_0_1_percent);

      InfoGEQFloatFilter info_1_filter(filter_field.value(), 0.01, false);
      filtered_population = filtered_population->filterVariants(info_1_filter);
      ExecEnv::log().info( "{}, Original VCF Population size: {}, Filtered Population size: {}",
                           info_1_filter.filterName(), vcf_population->variantCount(), filtered_population->variantCount());
      filtered_population->processAll(age_analysis_1_percent);

      InfoGEQFloatFilter info_2_filter(filter_field.value(), 0.02, false);
      filtered_population = filtered_population->filterVariants(info_2_filter);
      ExecEnv::log().info( "{}, Original VCF Population size: {}, Filtered Population size: {}",
                           info_2_filter.filterName(), vcf_population->variantCount(), filtered_population->variantCount());
      filtered_population->processAll(age_analysis_2_percent);

      InfoGEQFloatFilter info_5_filter(filter_field.value(), 0.05, false);
      filtered_population = filtered_population->filterVariants(info_5_filter);
      ExecEnv::log().info( "{}, Original VCF Population size: {}, Filtered Population size: {}",
                           info_5_filter.filterName(), vcf_population->variantCount(), filtered_population->variantCount());
      filtered_population->processAll(age_analysis_5_percent);

      InfoGEQFloatFilter info_10_filter(filter_field.value(), 0.10, false);
      filtered_population = filtered_population->filterVariants(info_10_filter);
      ExecEnv::log().info( "{}, Original VCF Population size: {}, Filtered Population size: {}",
                           info_10_filter.filterName(), vcf_population->variantCount(), filtered_population->variantCount());
      filtered_population->processAll(age_analysis_10_percent);

      ExecEnv::log().info("Write age statistics for vcf (population) file: {}", vcf_population->populationId());

      std::ofstream outfile(output_file_name_,  std::ofstream::out | std::ofstream::app);

      if (not outfile.good()) {

        ExecEnv::log().error("InfoFilterAnalysis::finalizeAnalysis; could not open results file: {}", output_file_name_);
        return false;

      }

      outfile << "Age analysis for vcf file (population), " << vcf_population->populationId() << '\n';
      outfile << "All Variants, " << age_analysis_all.variantCount() << '\n';
      outfile << age_analysis_all;
      outfile << "Filtered AF >= 0.1%, " << age_analysis_0_1_percent.variantCount() << '\n';
      outfile << age_analysis_0_1_percent;
      outfile << "Filtered AF >= 1%, " << age_analysis_1_percent.variantCount() << '\n';
      outfile << age_analysis_1_percent;
      outfile << "Filtered AF >= 2%, " << age_analysis_2_percent.variantCount() << '\n';
      outfile << age_analysis_2_percent;
      outfile << "Filtered AF >= 5%, " << age_analysis_5_percent.variantCount() << '\n';
      outfile << age_analysis_5_percent;
      outfile << "Filtered AF >= 10%, " << age_analysis_10_percent.variantCount() << '\n';
      outfile << age_analysis_10_percent << '\n' << '\n' << '\n';

      outfile.flush();

      age_analysis_all_.addAgeAnalysis(age_analysis_all);
      age_analysis_0_1_percent_.addAgeAnalysis(age_analysis_0_1_percent);
      age_analysis_1_percent_.addAgeAnalysis(age_analysis_1_percent);
      age_analysis_2_percent_.addAgeAnalysis(age_analysis_2_percent);
      age_analysis_5_percent_.addAgeAnalysis(age_analysis_5_percent);
      age_analysis_10_percent_.addAgeAnalysis(age_analysis_10_percent);


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

  std::ofstream outfile(output_file_name_,  std::ofstream::out | std::ofstream::app);

  if (not outfile.good()) {

    ExecEnv::log().error("InfoFilterAnalysis::finalizeAnalysis; could not open results file: {}", output_file_name_);
    return false;

  }

  outfile << "Age analysis for ALL vcf files (populations)" << '\n';
  outfile << "All Variants, " << age_analysis_all_.variantCount() << '\n';
  outfile << age_analysis_all_;
  outfile << "Filtered AF >= 0.1%, " << age_analysis_0_1_percent_.variantCount() << '\n';
  outfile << age_analysis_0_1_percent_;
  outfile << "Filtered AF >= 1%, " << age_analysis_1_percent_.variantCount() << '\n';
  outfile << age_analysis_1_percent_;
  outfile << "Filtered AF >= 2%, " << age_analysis_2_percent_.variantCount() << '\n';
  outfile << age_analysis_2_percent_;
  outfile << "Filtered AF >= 5%, " << age_analysis_5_percent_.variantCount() << '\n';
  outfile << age_analysis_5_percent_;
  outfile << "Filtered AF >= 10%, " << age_analysis_10_percent_.variantCount() << '\n';
  outfile << age_analysis_10_percent_;

  outfile.flush();

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
