//
// Created by kellerberrin on 20/11/20.
//

#ifndef KGL_ANALYSIS_INBREED_FREQ_H
#define KGL_ANALYSIS_INBREED_FREQ_H

#include "kgl_hsgenealogy_parser.h"
#include "kel_optimize.h"
#include "kgl_variant_db_freq.h"

namespace kellerberrin::genome {   //  organization::project level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold allele class frequencies based on the potential alleles occurring at a location.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AlleleClassFrequencies {

public:

  AlleleClassFrequencies(double major_hom, double major_het, double minor_hom, double minor_het, double inbreeding) :
      major_homozygous_(major_hom), major_heterozygous_(major_het), minor_homozygous_(minor_hom),
      minor_heterozygous_(minor_het), inbreeding_(inbreeding) {}

  ~AlleleClassFrequencies() = default;

  [[nodiscard]] bool validFrequencies() const { return std::fabs(sumFrequencies() - 1.0) < epsilon_; }

  [[nodiscard]] double sumFrequencies() const { return major_homozygous_ + major_heterozygous_ + minor_homozygous_ + minor_heterozygous_; }

  [[nodiscard]] double majorHomozygous() const { return major_homozygous_; }

  [[nodiscard]] double majorHeterozygous() const { return major_heterozygous_; }

  [[nodiscard]] double minorHomozygous() const { return minor_homozygous_; }

  [[nodiscard]] double minorHeterozygous() const { return minor_heterozygous_; }

  [[nodiscard]] double inbreeding() const { return inbreeding_; }

  void nonNegative() {

    major_homozygous_ = std::max(0.0, major_homozygous_);
    major_heterozygous_ = std::max(0.0, major_heterozygous_);
    minor_homozygous_ = std::max(0.0, minor_homozygous_);
    minor_heterozygous_ = std::max(0.0, minor_heterozygous_);

  }

  void normalize() {

    nonNegative();
    double sum_freqs = sumFrequencies();
    major_homozygous_ = major_homozygous_ / sum_freqs;
    major_heterozygous_ = major_heterozygous_ / sum_freqs;
    minor_homozygous_ = minor_homozygous_ / sum_freqs;
    minor_heterozygous_ = minor_heterozygous_ / sum_freqs;

  }


private:

  double major_homozygous_{0.0};
  double major_heterozygous_{0.0};
  double minor_homozygous_{0.0};
  double minor_heterozygous_{0.0};
  double inbreeding_{0.0};

  // Tolerance for the sum of allele classes
  constexpr static const double epsilon_{1.0e-04};

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MAJOR_HOMOZYGOUS - Two homozygous major alleles (both unrecorded).
// MAJOR_HETEROZYGOUS - One minor allele and the major allele (unrecorded).
// MINOR_HOMOZYGOUS - Two identical minor alleles.
// MINOR_HETEROZYGOUS - Two different heterozygous minor alleles (both recorded).
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enum class AlleleClassType {
  MAJOR_HOMOZYGOUS, MAJOR_HETEROZYGOUS, MINOR_HETEROZYGOUS, MINOR_HOMOZYGOUS
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures to hold allele frequency information at an offset.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AlleleFreqRecord {

public:

  AlleleFreqRecord(std::shared_ptr<const Variant> allele, double frequency)
      : allele_(std::move(allele)), frequency_(frequency) {}

  ~AlleleFreqRecord() = default;

  [[nodiscard]] const std::shared_ptr<const Variant> &allele() const { return allele_; }

  [[nodiscard]] double frequency() const { return frequency_; }

private:

  std::shared_ptr<const Variant> allele_;
  double frequency_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A vector of potential allele frequencies at a location.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AlleleFreqVector {

public:

  AlleleFreqVector( const OffsetDBArray& variant_vector,
                    const std::string& frequency_field);

  ~AlleleFreqVector() = default;

  // Return the allele vector.
  [[nodiscard]] const std::vector <AlleleFreqRecord> &alleleFrequencies() const { return allele_frequencies_; }

  // Sum of all minor alleles clamped to [0,1]
  [[nodiscard]] double minorAlleleFrequencies() const;

  // Complement of the sum of all minor alleles.
  [[nodiscard]] double majorAlleleFrequency() const;

  // Checks that we have a valid vector of minor alleles.
  // The minor alleles must sum to 1.0 or less and the allele vector must be non-empty.
  [[nodiscard]] bool checkValidAlleleVector();

  // Return the unadjusted allele class frequencies (may be -ve).
  [[nodiscard]] AlleleClassFrequencies unadjustedAlleleClassFrequencies(double inbreeding) const;

  // Return the allele class frequencies.
  [[nodiscard]] AlleleClassFrequencies alleleClassFrequencies(double inbreeding) const;

  // Randomly select an allele class outcome based on a unit [0, 1] random number.
  [[nodiscard]] AlleleClassType selectAlleleClass(double unit_rand, const AlleleClassFrequencies &class_freqs) const;

  // Randomly select a minor homozygous allele based on a unit [0,1] random number, std::nullopt if error (no minor allele).
  [[nodiscard]] std::optional <AlleleFreqRecord> selectMinorHomozygous(double unit_rand, const AlleleClassFrequencies &class_freqs) const;

  // Randomly select a minor heterozygous allele based on a unit [0,1] random number, std::nullopt if error (no minor allele).
  [[nodiscard]] std::optional <AlleleFreqRecord> selectMajorHeterozygous(double unit_rand, const AlleleClassFrequencies &class_freqs) const;

  // Randomly select a pair of distinct minor alleles based on two random numbers, std::nullopt if error (not two minor alleles).
  [[nodiscard]] std::optional <std::pair<AlleleFreqRecord, AlleleFreqRecord>>
  selectMinorHeterozygous(double unit_rand, const AlleleClassFrequencies &class_freqs) const;

private:

  // All minor allele frequencies at minor allele offset (from Gnomad)
  std::vector <AlleleFreqRecord> allele_frequencies_;
  // Tolerance for the sum of alleles
  constexpr static const double epsilon_class_{1.0e-5};

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
                 const AlleleFreqRecord &first_allele_freq,
                 const AlleleFreqRecord &second_allele_freq,
                 const AlleleFreqVector &allele_frequencies) :
      allele_type_(allele_type),
      first_allele_freq_(first_allele_freq),
      second_allele_freq_(second_allele_freq),
      allele_frequencies_(allele_frequencies) {}

  ~AlleleFreqInfo() = default;

  [[nodiscard]] AlleleClassType alleleType() const { return allele_type_; }

  // The probfailure of the allele type specified above.
  [[nodiscard]] const AlleleFreqRecord &firstAllele() const { return first_allele_freq_; }

  [[nodiscard]] const AlleleFreqRecord &secondAllele() const { return second_allele_freq_; }

  [[nodiscard]] const AlleleFreqVector &alleleFrequencies() const { return allele_frequencies_; }

private:

  // Actual occurrence of the allele type and frequencies.
  const AlleleClassType allele_type_;
  AlleleFreqRecord first_allele_freq_;
  AlleleFreqRecord second_allele_freq_;
  // All the minor alleles defined for this location.
  // From the reference variant genome.
  AlleleFreqVector allele_frequencies_;

};


} // namespace.




#endif //KGL_ANALYSIS_MUTATION_INBREED_FREQ_H
