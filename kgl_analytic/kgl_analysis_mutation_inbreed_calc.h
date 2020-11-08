//
// Created by kellerberrin on 27/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_CALC_H
#define KGL_ANALYSIS_MUTATION_INBREED_CALC_H


#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"
#include "kel_optimize.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// Called by the threadpool for each genome/sample.
struct LocusResults {

  GenomeId_t genome;
  size_t major_hetero_count{0};   // 1 minor allele only.
  double major_hetero_freq{0.0};
  size_t minor_hetero_count{0};   // 2 different minor alleles.
  double minor_hetero_freq{0.0};
  size_t minor_homo_count{0};     // 2 identical minor alleles.
  double minor_homo_freq{0.0};
  size_t major_homo_count{0};     // 2 identical major alleles (generally not recorded).
  double major_homo_freq{0.0};
  size_t total_allele_count{0};  // All alleles.
  double inbred_allele_sum{0.0};

};

// MINOR_HOMOZYGOUS - Two identical minor alleles.
// MAJOR_HETEROZYGOUS - One minor allele and the major allele (unrecorded).
// MINOR_HETEROZYGOUS - Two different heterozygous minor alleles (both recorded).
// MAJOR_HOMOZYGOUS - Two homozygous major alleles (both unrecorded).
enum class MinorAlleleType { MINOR_HOMOZYGOUS, MAJOR_HETEROZYGOUS, MINOR_HETEROZYGOUS, MAJOR_HOMOZYGOUS };
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold allele frequency information at an offset.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AlleleFreqRecord {

public:

  AlleleFreqRecord(std::shared_ptr<const Variant> allele, double frequency) : allele_(std::move(allele)), frequency_(frequency) {}
  ~AlleleFreqRecord() = default;

  [[nodiscard]] const std::shared_ptr<const Variant>& allele() const { return allele_; }
  [[nodiscard]] double frequency() const { return frequency_; }

private:

  std::shared_ptr<const Variant> allele_;
  const double frequency_;

};

class AlleleFreqVector {

public:

  explicit AlleleFreqVector(std::vector<AlleleFreqRecord> allele_frequencies) : allele_frequencies_(std::move(allele_frequencies)) {}
  ~AlleleFreqVector() = default;

  [[nodiscard]] const std::vector<AlleleFreqRecord>& alleleFrequencies() const { return allele_frequencies_; }
  // Complement of the sum of all minor alleles.
  [[nodiscard]] double majorAlleleFrequency() const;
  // The probability of major homozygous alleles (unrecorded)
  [[nodiscard]] double majorHomozygous() const;
  // The probability of major heterozygous alleles, this is a single minor allele and major allele (unrecorded)
  [[nodiscard]] double majorHeterozygous() const;
  // The probability of minor allele homozygous alleles
  [[nodiscard]] double minorHomozygous() const;
  // The probability of minor and minor allele heterozygous alleles
  [[nodiscard]] double minorHeterozygous() const;
  // Sum of the allele frequencies.
  [[nodiscard]] double sumFrequencies() const;
  // True if duplicate alleles found
  [[nodiscard]] bool checkDuplicates() const;

private:

  // All minor allele frequencies at minor allele offset (from Gnomad)
  std::vector<AlleleFreqRecord> allele_frequencies_;

};



class AlleleFreqInfo {
public:

  AlleleFreqInfo(MinorAlleleType allele_type,
                 double first_allele_freq,
                 double second_allele_freq,
                 const AlleleFreqVector& allele_frequencies) :
      allele_type_(allele_type),
      first_allele_freq_(first_allele_freq),
      second_allele_freq_(second_allele_freq),
      allele_frequencies_(allele_frequencies) {}
  ~AlleleFreqInfo() = default;

  [[nodiscard]] MinorAlleleType alleleType() const { return allele_type_; }
  [[nodiscard]] double firstAlleleFrequency() const { return first_allele_freq_; }
  [[nodiscard]] double secondAlleleFrequency() const { return second_allele_freq_; }
  [[nodiscard]] const AlleleFreqVector& alleleFrequencies() const { return allele_frequencies_; }

private:

  const MinorAlleleType allele_type_;
  double first_allele_freq_{0.0};
  double second_allele_freq_{0.0};
  // From the reference variant genome (generally Gnomad)
  // All minor allele frequencies at minor allele offset (from Gnomad)
  AlleleFreqVector allele_frequencies_;

};

// Algorithm type to pass to the threadpool.
using InbreedingAlgorithm = std::function<LocusResults(const GenomeId_t& genome_id,
                                                       const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                       const std::string& super_population_field,
                                                       const std::shared_ptr<const ContigVariant>& locus_list)>;

// Static class to hold the inbreeding algorithms.
class InbreedingCalculation  {

public:

  InbreedingCalculation() = delete;
  ~InbreedingCalculation() = delete;


  // The inbreeding Algorithms
  [[nodiscard]] static LocusResults processRitlandLocus( const ContigId_t& genome_id,
                                                         const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                         const std::string& super_population_field,
                                                         const std::shared_ptr<const ContigVariant>& locus_list);
  
  [[nodiscard]] static LocusResults processSimple(const ContigId_t& contig_id,
                                                  const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                  const std::string& super_population_field,
                                                  const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static LocusResults processHallME( const GenomeId_t& genome_id,
                                                  const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                  const std::string& super_population_field,
                                                  const std::shared_ptr<const ContigVariant>& locus_list);

  // The loglikelihood algorithm.
  [[nodiscard]] static LocusResults processLogLikelihood(const GenomeId_t& genome_id,
                                                         const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                         const std::string& super_population_field,
                                                         const std::shared_ptr<const ContigVariant>& locus_list);

  // Algorithm key used in the the algorithm selection map.
  inline static const std::string RITLAND_LOCUS{"RitlandLocus"};
  inline static const std::string SIMPLE{"Simple"};
  inline static const std::string HALL_ME{"HallME"};
  inline static const std::string LOGLIKELIHOOD{"Loglikelihood"};

  static const std::map<std::string, InbreedingAlgorithm>& algoMap() { return inbreeding_algo_map_; }

private:

  inline static const std::map<std::string, InbreedingAlgorithm> inbreeding_algo_map_ = {
      {RITLAND_LOCUS, processRitlandLocus},
      {SIMPLE, processSimple},
      {HALL_ME, processHallME},
      {LOGLIKELIHOOD, processLogLikelihood}
  };

  // The initial guess for the Hall expectation maximization algorithm.
  constexpr static const double FINAL_ACCURACY_ = 1E-04;
  constexpr static const size_t MINIMUM_ITERATIONS_ = 50;
  constexpr static const size_t MAXIMUM_ITERATIONS_ = 1000;
  constexpr static const size_t MAX_RETRIES_ = 50;
  constexpr static const size_t MIN_RETRIES_ = 2;

  [[nodiscard]] static Optimize createLogLikelihoodOptimizer();

  [[nodiscard]] static std::pair<std::vector<AlleleFreqInfo>, LocusResults>
  generateFrequencies(const GenomeId_t& genome_id,
                      const std::shared_ptr<const DiploidContig>& contig_ptr,
                      const std::string& super_population_field,
                      const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static double logLikelihood(std::vector<double>& x, std::vector<AlleleFreqInfo>& data);

};



} // namespace















#endif //KGL_ANALYSIS_MUTATION_INBREED_CALC_H
