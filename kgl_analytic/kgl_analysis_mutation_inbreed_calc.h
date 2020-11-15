//
// Created by kellerberrin on 27/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_CALC_H
#define KGL_ANALYSIS_MUTATION_INBREED_CALC_H


#include "kgl_variant_db_phased.h"
#include "kgl_ped_parser.h"
#include "kel_optimize.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold summary actual and implied allele frequency information at an offset.
// Called by the threadpool for each genome/sample.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold allele class frequencies based on the potential alleles occurring at a location.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AlleleClassFrequencies {

public:

  AlleleClassFrequencies(double major_hom, double major_het, double minor_hom, double minor_het) :
      major_homozygous_(major_hom), major_heterozygous_(major_het), minor_homozygous_(minor_hom), minor_heterozygous_(minor_het) {}
  ~AlleleClassFrequencies() = default;

  [[nodiscard]] bool validFrequencies() const { return std::fabs(sumFrequencies() - 1.0) < epsilon_; }
  [[nodiscard]] double sumFrequencies() const { return major_homozygous_ + major_heterozygous_ + minor_homozygous_ + minor_heterozygous_; }
  [[nodiscard]] double majorHomozygous() const { return major_homozygous_; }
  [[nodiscard]] double majorHeterozygous() const { return major_heterozygous_; }
  [[nodiscard]] double minorHomozygous() const { return minor_homozygous_; }
  [[nodiscard]] double minorHeterozygous() const { return minor_heterozygous_; }

private:

  double major_homozygous_{0.0};
  double major_heterozygous_{0.0};
  double minor_homozygous_{0.0};
  double minor_heterozygous_{0.0};

  // Tolerance for the sum of allele classes
  constexpr static const double epsilon_{1.0e-04};

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAJOR_HOMOZYGOUS - Two homozygous major alleles (both unrecorded).
// MAJOR_HETEROZYGOUS - One minor allele and the major allele (unrecorded).
// MINOR_HOMOZYGOUS - Two identical minor alleles.
// MINOR_HETEROZYGOUS - Two different heterozygous minor alleles (both recorded).
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enum class AlleleClassType { MAJOR_HOMOZYGOUS, MAJOR_HETEROZYGOUS, MINOR_HETEROZYGOUS, MINOR_HOMOZYGOUS };

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold allele frequency information at an offset.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AlleleFreqRecord {

public:

  AlleleFreqRecord(std::shared_ptr<const Variant> allele, double frequency, std::string freq_field)
    : allele_(std::move(allele)), frequency_(frequency), freq_field_(std::move(freq_field)) {}
  ~AlleleFreqRecord() = default;

  [[nodiscard]] const std::shared_ptr<const Variant>& allele() const { return allele_; }
  [[nodiscard]] double frequency() const { return frequency_; }
  [[nodiscard]] const std::string& freqField() const { return freq_field_; }
  void frequency(double freq) { frequency_ = freq; }

private:

  std::shared_ptr<const Variant> allele_;
  double frequency_;
  std::string freq_field_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A vector of potential allele frequencies at a location.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AlleleFreqVector {

public:

  AlleleFreqVector(const std::vector<std::shared_ptr<const Variant>>& variant_vector, const std::string& frequency_field);
  ~AlleleFreqVector() = default;

  // Return the allele vector.
  [[nodiscard]] const std::vector<AlleleFreqRecord>& alleleFrequencies() const { return allele_frequencies_; }
  // Sum of all minor alleles.
  [[nodiscard]] double minorAlleleFrequencies() const;
  // Complement of the sum of all minor alleles.
  [[nodiscard]] double majorAlleleFrequency() const;
  // Checks that we have a valid vector of minor alleles.
  // Updates allele frequencies under limited circumstances.
  [[nodiscard]] bool checkValidAlleleVector();
  // Return the allele class frequencies.
  [[nodiscard]] AlleleClassFrequencies alleleClassFrequencies(double inbreeding) const;
  // Randomly select an allele class outcome based on a unit [0, 1] random number.
  [[nodiscard]] AlleleClassType selectAlleleClass(double unit_rand, double inbreeding) const;
  // Randomly select a minor allele based on a unit [0,1] random number, std::nullopt if error (no minor allele).
  [[nodiscard]] std::optional<AlleleFreqRecord> selectMinorAllele(double unit_rand) const;
  // Randomly select a pair of distinct minor alleles based on two random numbers, std::nullopt if error (not two minor alleles).
  [[nodiscard]] std::optional<std::pair<AlleleFreqRecord, AlleleFreqRecord>> selectPairMinorAllele(double unit_rand1, double unit_rand2) const;
  // The probability of major homozygous alleles (unrecorded)
  [[nodiscard]] double majorHomozygous(double inbreeding) const;
  // The probability of major heterozygous alleles, this is a single minor allele and major allele (unrecorded)
  [[nodiscard]] double majorHeterozygous(double inbreeding) const;
  // The probability of minor allele homozygous alleles
  [[nodiscard]] double minorHomozygous(double inbreeding) const;
  // The probability of minor and minor allele heterozygous alleles
  [[nodiscard]] double minorHeterozygous(double inbreeding) const;

private:

  // All minor allele frequencies at minor allele offset (from Gnomad)
  std::vector<AlleleFreqRecord> allele_frequencies_;
  // Tolerance for the sum of allele classes
  constexpr static const double epsilon_class_{1.0e-10};
  // Tolerance for the sum of allele frequencies
  constexpr static const double epsilon_sum_{1.0e-03};

  // True if duplicate alleles found
  [[nodiscard]] bool checkDuplicates() const;
  // Not clamped to 0.0, 1.0
  [[nodiscard]] double sumAlleleFrequencies() const;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The actual allele class at a location, the actual two allele frequencies, and a vector of potential allele frequencies at a location.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AlleleFreqInfo {

public:

  AlleleFreqInfo(AlleleClassType allele_type,
                 const AlleleFreqRecord& first_allele_freq,
                 const AlleleFreqRecord& second_allele_freq,
                 const AlleleFreqVector& allele_frequencies) :
      allele_type_(allele_type),
      first_allele_freq_(first_allele_freq),
      second_allele_freq_(second_allele_freq),
      allele_frequencies_(allele_frequencies) {}
  ~AlleleFreqInfo() = default;

  [[nodiscard]] AlleleClassType alleleType() const { return allele_type_; }
  // The probability of the allele type specified above.
  [[nodiscard]] double alleleTypeFrequency(double inbreeding) const;
  [[nodiscard]] const AlleleFreqRecord& firstAllele() const { return first_allele_freq_; }
  [[nodiscard]] const AlleleFreqRecord& secondAllele() const { return second_allele_freq_; }
  [[nodiscard]] const AlleleFreqVector& alleleFrequencies() const { return allele_frequencies_; }

private:

  // Actual occurrence of the allele type and frequencies.
  const AlleleClassType allele_type_;
  AlleleFreqRecord first_allele_freq_;
  AlleleFreqRecord second_allele_freq_;
  // All the minor alleles defined for this location.
  // From the reference variant genome.
  AlleleFreqVector allele_frequencies_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////
// The Functional definition of an inbreeding algorithm pass to the threadpool
////////////////////////////////////////////////////////////////////////////////////////////////

// Algorithm type to pass to the threadpool.
using InbreedingAlgorithm = std::function<LocusResults(const GenomeId_t& genome_id,
                                                       const std::shared_ptr<const DiploidContig>& contig_ptr,
                                                       const std::string& super_population_field,
                                                       const std::shared_ptr<const ContigVariant>& locus_list)>;

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
//      {SIMPLE, processSimple},
//      {HALL_ME, processHallME},
//      {LOGLIKELIHOOD, processLogLikelihood}
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
