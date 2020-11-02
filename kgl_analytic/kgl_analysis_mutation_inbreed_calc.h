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
  size_t major_hetero_count{0};
  size_t minor_hetero_count{0};
  size_t homo_count{0};
  size_t total_allele_count{0};
  double inbred_allele_sum{0.0};

};

// HOMOZYGOUS - Two identical minor alleles.
// HETEROZYGOUS - One minor allele and the major allele (unrecorded).
// MINOR_HETEROZYGOUS - Two different heterozygous minor alleles (both recorded).
enum class MinorAlleleType { HOMOZYGOUS, HETEROZYGOUS, MINOR_HETEROZYGOUS };
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


class AlleleFreqInfo {
public:

  AlleleFreqInfo(MinorAlleleType allele_type,
                 const AlleleFreqRecord& minor_allele,
                 const AlleleFreqRecord& second_allele,
                 std::vector<AlleleFreqRecord>&& allele_frequencies) :
      allele_type_(allele_type),
      minor_allele_(minor_allele),
      second_allele_(second_allele),
      allele_frequencies_(allele_frequencies) {}
  ~AlleleFreqInfo() = default;

  [[nodiscard]] MinorAlleleType alleleType() const { return allele_type_; }
  [[nodiscard]] const AlleleFreqRecord& minorAllele() const { return minor_allele_; }
  [[nodiscard]] const AlleleFreqRecord& secondAllele() const { return second_allele_; }
  [[nodiscard]] const std::vector<AlleleFreqRecord>& alleleFrequencies() const { return allele_frequencies_; }
  [[nodiscard]] double majorAlleleFrequency() const;

private:

  const MinorAlleleType allele_type_;
  // From the diploid genome.
  const AlleleFreqRecord minor_allele_;
  // From the diploid genome.
  // If homozygous then minor allele, else major allele freq (same variant_ptr as minor), else different minor allele (rare).
  // Warning, do not compare with minor allele to determine if homozygous, as the major allele is represented by the
  // minor allele variant pointer. (Note to self, do we need a major allele reference = alternate allele variant type?)
  const AlleleFreqRecord second_allele_;
  // From the reference variant genome (generally Gnomad)
  // All minor allele frequencies at minor allele offset (from Gnomad)
  std::vector<AlleleFreqRecord> allele_frequencies_;

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
  
  [[nodiscard]] static LocusResults processRitlandMME( const ContigId_t& contig_id,
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
  inline static const std::string RITLAND_MME{"RitlandMME"};
  inline static const std::string HALL_ME{"HallME"};
  inline static const std::string LOGLIKELIHOOD{"Loglikelihood"};

  static const std::map<std::string, InbreedingAlgorithm>& algoMap() { return inbreeding_algo_map_; }

private:

  inline static const std::map<std::string, InbreedingAlgorithm> inbreeding_algo_map_ = {
      {RITLAND_LOCUS, processRitlandLocus},
      {RITLAND_MME, processRitlandMME},
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
  generateGnomadFreq(const GenomeId_t& genome_id,
                     const std::shared_ptr<const DiploidContig>& contig_ptr,
                     const std::string& super_population_field,
                     const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static std::pair<std::vector<AlleleFreqInfo>, LocusResults>
  generateDiploidFreq(const GenomeId_t& genome_id,
                     const std::shared_ptr<const DiploidContig>& contig_ptr,
                     const std::string& super_population_field,
                     const std::shared_ptr<const ContigVariant>& locus_list);

  [[nodiscard]] static double logLikelihood(std::vector<double>& x, std::vector<AlleleFreqInfo>& data);

};



} // namespace















#endif //KGL_ANALYSIS_MUTATION_INBREED_CALC_H
