//
// Created by kellerberrin on 27/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_CALC_H
#define KGL_ANALYSIS_MUTATION_INBREED_CALC_H


#include "kgl_analysis_mutation_inbreed_freq.h"
#include "kgl_analysis_mutation_inbreed_locus.h"

#include <list>

namespace kellerberrin::genome {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parameter object passed to the inbreeding analysis algorithms.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class InbreedingParameters {

public:

  InbreedingParameters() = default;

  ~InbreedingParameters() = default;

  [[nodiscard]] const LociiVectorArguments &lociiArguments() const { return locii_selection_; }

  // Non const version.
  [[nodiscard]] LociiVectorArguments &lociiArguments() { return locii_selection_; }

  [[nodiscard]] const std::string &inbreedingAlgorthim() const { return inbreeding_algorithm_; }

  [[nodiscard]] const std::string &outputFile() const { return output_file_; }

  void lociiArguments(const LociiVectorArguments &locii_args) { locii_selection_ = locii_args; }

  void inbreedingAlgorthim(const std::string &algo_name) { inbreeding_algorithm_ = algo_name; }

  void outputFile(const std::string &output_file) { output_file_ = output_file; }

private:

  LociiVectorArguments locii_selection_;
  std::string inbreeding_algorithm_{"Loglikelihood"};
  VariantDatabaseSource variant_source_{VariantDatabaseSource::GNOMAD2_1};
  std::string output_file_{"inbreeding_output"};

};


/////////////////////////////////////////////////////////////////////////////////////////////////
// The Functional definition of an inbreeding algorithm pass to the threadpool
////////////////////////////////////////////////////////////////////////////////////////////////

// Algorithm type to pass to the threadpool.
using InbreedingAlgorithm = std::function<LocusResults(const GenomeId_t& genome_id,
                                                       const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                       const std::string& super_population_field,
                                                       const std::shared_ptr<const ContigVariant>& locus_list,
                                                       const InbreedingParameters& parameters)>;

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple class the check retries
//
///////////////////////////////////////////////////////////////////////////////////////////////

class RetryCalcResult {

public:

  RetryCalcResult(double tolerance, size_t min_retry, size_t max_retry) : tolerance_(tolerance), min_retry_(min_retry), max_retry_(max_retry) {}
  ~RetryCalcResult() = default;

  // Returns true if retries >= min_try and tolerance OK or if max_retry is exceeded.
  bool checkRetry(double retry);
  size_t retries() const { return retry_count_; }

private:

  const double tolerance_;
  const size_t min_retry_;
  const size_t max_retry_;
  size_t retry_count_{0};
  std::list<double> current_retries_;

  bool checkTolerance() const;

};



////////////////////////////////////////////////////////////////////////////////////////////////
// Static class to hold the inbreeding algorithms.
////////////////////////////////////////////////////////////////////////////////////////////////

class InbreedingCalculation  {

public:

  InbreedingCalculation() = delete;
  ~InbreedingCalculation() = delete;


  // The inbreeding Algorithms
  [[nodiscard]] static LocusResults processRitlandLocus( const ContigId_t& genome_id,
                                                         const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                         const std::string& super_population_field,
                                                         const std::shared_ptr<const ContigVariant>& locus_list,
                                                         const InbreedingParameters& parameters);
  
  [[nodiscard]] static LocusResults processSimple(const ContigId_t& contig_id,
                                                  const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                  const std::string& super_population_field,
                                                  const std::shared_ptr<const ContigVariant>& locus_list,
                                                  const InbreedingParameters& parameters);

  [[nodiscard]] static LocusResults processHallME( const GenomeId_t& genome_id,
                                                  const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                  const std::string& super_population_field,
                                                  const std::shared_ptr<const ContigVariant>& locus_list,
                                                  const InbreedingParameters& parameters);

  // The loglikelihood algorithm.
  [[nodiscard]] static LocusResults processLogLikelihood(const GenomeId_t& genome_id,
                                                         const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                         const std::string& super_population_field,
                                                         const std::shared_ptr<const ContigVariant>& locus_list,
                                                         const InbreedingParameters& parameters);

  // Algorithm key used in the the algorithm selection map.
  inline static const std::string RITLAND_LOCUS{"RitlandLocus"};
  inline static const std::string SIMPLE{"Simple"};
  inline static const std::string HALL_ME{"HallME"};
  inline static const std::string LOGLIKELIHOOD{"Loglikelihood"};

  static const std::map<std::string, InbreedingAlgorithm>& algoMap() { return inbreeding_algo_map_; }
  static std::optional<InbreedingAlgorithm> namedAlgorithm(const std::string& algorithm_name);

private:

  inline static const std::map<std::string, InbreedingAlgorithm> inbreeding_algo_map_ = {
      {RITLAND_LOCUS, processRitlandLocus},
      {SIMPLE, processSimple},
      {HALL_ME, processHallME},
      {LOGLIKELIHOOD, processLogLikelihood}
  };

  // The initial guess for the Hall expectation maximization algorithm.
  constexpr static const double FINAL_ACCURACY_ = 1E-04;
  constexpr static const double INIT_UPPER_ = 0.5;
  constexpr static const double INIT_LOWER_ = -0.5;
  constexpr static const size_t MINIMUM_ITERATIONS_ = 50;
  constexpr static const size_t MAXIMUM_ITERATIONS_ = 1000;
  constexpr static const size_t MAX_RETRIES_ = 50;
  constexpr static const size_t MIN_RETRIES_ = 5;

  [[nodiscard]] static Optimize createLogLikelihoodOptimizer();

  [[nodiscard]] static std::pair<std::vector<AlleleFreqInfo>, LocusResults>
  generateFrequencies(const GenomeId_t& genome_id,
                      const std::shared_ptr<const DiploidContig>& contig_ptr,
                      const std::string& super_population_field,
                      const std::shared_ptr<const ContigVariant>& locus_list,
                      VariantDatabaseSource variantSource);

  [[nodiscard]] static double logLikelihood(std::vector<double>& x, std::vector<AlleleFreqInfo>& data);

};



} // namespace















#endif //KGL_ANALYSIS_MUTATION_INBREED_CALC_H
