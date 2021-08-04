//
// Created by kellerberrin on 30/12/20.
//

#ifndef KGL_ANALYSIS_INBREED_ARGS_H
#define KGL_ANALYSIS_INBREED_ARGS_H


#include "kel_exec_env.h"
#include "kgl_data_file_type.h"
#include "kgl_runtime.h"

#include <string>

namespace kellerberrin::genome {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple convenience class to decode XML parsed arguments into a suitable data structure.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The current argument XML structure is a vector of the following.
//
//<help>Generate the inbreeding coefficient using "Synthetic" or real "Inbreed" data</help>
//<parameterIdent>AnalysisType</parameterIdent>
//
//<help>Output file name; Text</help>
//<parameterIdent>OutputFile</parameterIdent>
//
//<help>Minimum Allele Frequency to use in the analysis; Float</help>
//<parameterIdent>MinAlleleFreq</parameterIdent>
//
//<help>Maximum Allele Frequency to use in the analysis; Float</help>
//<parameterIdent>MaxAlleleFreq</parameterIdent>
//
//<help>The Inbreeding Algorithm to Use; Text</help>
//<parameterIdent>Algorithm</parameterIdent>
//
//<help>The Minimum Contig Offset to Examine; Integer</help>
//<parameterIdent>LowerWindow</parameterIdent>
//
//<help>The Maximum Contig Offset to examine, set to High Value for entire contig; Integer</help>
//<parameterIdent>UpperWindow</parameterIdent>
//
//<help>The Number of allele offsets to examine, set to High Value for entire contig; Integer</help>
//<parameterIdent>LociiCount</parameterIdent>
//
//<help>The minimum distance between allele offsets to minimize any linkage disequilibrium effects; Integer</help>
//<parameterIdent>SamplingDistance</parameterIdent>
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passes arguments remotely to the locii vector class
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class LociiVectorArguments {

public:


  LociiVectorArguments() =default;
  ~LociiVectorArguments() = default;

  [[nodiscard]] ContigOffset_t lowerOffset() const { return lower_offset_; }
  [[nodiscard]] ContigOffset_t upperOffset() const { return upper_offset_; }
  [[nodiscard]] size_t lociiSpacing() const { return spacing_; };
  [[nodiscard]] size_t lociiCount () const { return locii_count; }
  [[nodiscard]] double minAlleleFrequency() const { return allele_frequency_min_; }
  [[nodiscard]] double maxAlleleFrequency() const { return allele_frequency_max_; }
  [[nodiscard]] DataSourceEnum frequencySource() const { return frequency_source_; }

  void lowerOffset(ContigOffset_t lower) { lower_offset_ = lower; }
  void upperOffset(ContigOffset_t upper) { upper_offset_ = upper; }
  void lociiSpacing(size_t spacing) { spacing_ = spacing; };
  void lociiCount (size_t count) { locii_count = count; }
  void minAlleleFrequency(double min_AF) { allele_frequency_min_ = std::clamp(min_AF, 0.0, 1.0); }
  void maxAlleleFrequency(double max_AF) { allele_frequency_max_ = std::clamp(max_AF, 0.0, 1.0); }
  void frequencySource(DataSourceEnum frequency_source) { frequency_source_ = frequency_source; }

private:

  ContigOffset_t lower_offset_{0};
  ContigOffset_t upper_offset_{1000000000};
  size_t spacing_{1000};
  size_t locii_count{1000};
  double allele_frequency_min_{0.0};
  double allele_frequency_max_{1.0};
  DataSourceEnum frequency_source_{DataSourceEnum::Gnomad2_1};

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passes arguments remotely to the execute class
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class InbreedingParameters {

public:

  InbreedingParameters() = default;
  ~InbreedingParameters() = default;

  const std::string& paramIdent() const { return parameter_ident_; }
  void paramIdent(const std::string& ident) { parameter_ident_ = ident; }

  [[nodiscard]] const LociiVectorArguments &lociiArguments() const { return locii_selection_; }
  // Non const version.
  [[nodiscard]] LociiVectorArguments &lociiArguments() { return locii_selection_; }

  [[nodiscard]] const std::string &inbreedingAlgorthim() const { return inbreeding_algorithm_; }
  void inbreedingAlgorthim(const std::string &algo_name) { inbreeding_algorithm_ = algo_name; }

  // Create synthetic data to analyze or analyze actual supplied genome data.
  [[nodiscard]] bool analyzeSynthetic() const { return analyze_synthetic_; }
  void analyzeSynthetic(bool synthetic) {  analyze_synthetic_ = synthetic; }

  // Output file name.
  [[nodiscard]] const std::string& outputFile() const { return output_file_; }
  void outputFile(const std::string& file_name) {  output_file_ = file_name; }


private:

  LociiVectorArguments locii_selection_;
  std::string parameter_ident_{"ParamIdent"};
  std::string inbreeding_algorithm_{"Loglikelihood"};
  std::string output_file_{"output"};
  bool analyze_synthetic_{true};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Extracts XML argument format into a vector of InbreedingParameters
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class InbreedArguments {

public:

  InbreedArguments() = delete;
  ~InbreedArguments() = delete;

  static std::vector<InbreedingParameters> extractParameters(const ActiveParameterList& named_parameters);

private:

  inline static const std::string ANALYSISTYPE_{"AnalysisType"};
  inline static const std::pair<std::string, size_t> OUTPUTFILE_{"OutputFile", 1};
  inline static const std::pair<std::string, size_t> MINALLELEFREQ_{"MinAlleleFreq", 1};
  inline static const std::pair<std::string, size_t> MAXALLELEFREQ_{"MaxAlleleFreq", 1};
  inline static const std::pair<std::string, size_t> ALGORITHM_{"Algorithm", 1};
  inline static const std::pair<std::string, size_t> LOWERWINDOW_{"LowerWindow", 1};
  inline static const std::pair<std::string, size_t> UPPERWINDOW_{"UpperWindow", 1};
  inline static const std::pair<std::string, size_t> LOCIICOUNT_{"LociiCount", 1};
  inline static const std::pair<std::string, size_t> SAMPLINGDISTANCE_{"SamplingDistance", 1};


};


} // namespace


#endif //KGL_ANALYSIS_INBREED_ARGS_H


