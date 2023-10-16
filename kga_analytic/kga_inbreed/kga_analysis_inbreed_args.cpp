//
// Created by kellerberrin on 30/12/20.
//

#include "kga_analysis_inbreed_args.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;

// This function extracts the values from the parsed XML file and constructs the analysis arguments.
std::vector<kgl::InbreedingParameters> kgl::InbreedArguments::extractParameters(const ActiveParameterList& named_parameters) {


  std::vector<kgl::InbreedingParameters> param_vector;

  for (auto const& named_block : named_parameters.getMap()) {

    auto [block_name, block_vector] = named_block.second;

    for (auto const& xml_vector : block_vector) {

      InbreedingParameters parameters;
      parameters.paramIdent(block_name);

      auto analysis_opt = xml_vector.getBool(ANALYSISTYPE_);
      if (analysis_opt) {

        bool synthetic_flag = analysis_opt.value();
        parameters.analyzeSynthetic(synthetic_flag);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for synthetic flag: {}", ANALYSISTYPE_);
        continue;

      }

      auto output_opt = xml_vector.getString(OUTPUTFILE_);
      if (output_opt) {

        std::string output_file = output_opt.value().front();
        parameters.outputFile(output_file);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for synthetic flag: {}", OUTPUTFILE_.first);
        continue;

      }

      auto algo_opt = xml_vector.getString(ALGORITHM_);
      if (algo_opt) {

        std::string algorithm = algo_opt.value().front();
        parameters.inbreedingAlgorthim(algorithm);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for inbreeding algorithm: {}", ALGORITHM_.first);
        continue;

      }

      auto min_freq_opt = xml_vector.getFloat(MINALLELEFREQ_);
      if (min_freq_opt) {

        double min_freq = min_freq_opt.value().front();
        parameters.lociiArguments().minAlleleFrequency(min_freq);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for minimum allele frequency: {}", MINALLELEFREQ_.first);
        continue;

      }

      auto max_freq_opt = xml_vector.getFloat(MAXALLELEFREQ_);
      if (max_freq_opt) {

        double max_freq = max_freq_opt.value().front();
        parameters.lociiArguments().maxAlleleFrequency(max_freq);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for minimum allele frequency: {}", MAXALLELEFREQ_.first);
        continue;

      }

      auto low_window_opt = xml_vector.getSize(LOWERWINDOW_);
      if (low_window_opt) {

        size_t lower_window = low_window_opt.value().front();
        parameters.lociiArguments().lowerOffset(lower_window);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for lower window offset: {}", LOWERWINDOW_.first);
        continue;

      }

      auto high_window_opt = xml_vector.getSize(UPPERWINDOW_);
      if (high_window_opt) {

        size_t upper_window = high_window_opt.value().front();
        parameters.lociiArguments().upperOffset(upper_window);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for upper window offset: {}", UPPERWINDOW_.first);
        continue;

      }

      auto locii_opt = xml_vector.getSize(LOCIICOUNT_);
      if (locii_opt) {

        size_t locii_count = locii_opt.value().front();
        parameters.lociiArguments().lociiCount(locii_count);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for locii count: {}", LOCIICOUNT_.first);
        continue;

      }

      auto sampling_opt = xml_vector.getSize(SAMPLINGDISTANCE_);
      if (sampling_opt) {

        size_t sampling_distance = sampling_opt.value().front();
        parameters.lociiArguments().lociiSpacing(sampling_distance);

      } else {

        ExecEnv::log().error("InbreedArguments::extractParameters; bad value for locii count: {}", SAMPLINGDISTANCE_.first);
        continue;

      }

      param_vector.push_back(parameters);

    }

    ExecEnv::log().info("Inbreeding Analysis, parsed named argument vector: {}, size: {}", block_name, param_vector.size());

  }

  return param_vector;

}
