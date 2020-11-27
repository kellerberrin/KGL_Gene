//
// Created by kellerberrin on 24/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_AUX_H
#define KGL_ANALYSIS_MUTATION_INBREED_AUX_H

#include "kgl_variant_db_phased.h"
#include "kgl_analysis_mutation_inbreed_calc.h"
#include "kgl_analysis_mutation_inbreed_freqdb.h"

#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passes arguments remotely to the locii vector class
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class LociiType { LOCII_FROM_TO, LOCII_FROM_COUNT };

class LociiVectorArguments {

public:


  explicit LociiVectorArguments(std::shared_ptr<const ContigVariant> unphased_contig_ptr) : unphased_contig_ptr_(unphased_contig_ptr) {}
  ~LociiVectorArguments() = default;

  [[nodiscard]] LociiType lociiType() const { return locii_type_; }
  [[nodiscard]] ContigOffset_t lowerOffset() const { return lower_offset_; }
  [[nodiscard]] ContigOffset_t upperOffset() const { return upper_offset_; }
  [[nodiscard]] size_t lociiSpacing() const { return spacing_; };
  [[nodiscard]] size_t lociiCount () const { return locii_count; }
  [[nodiscard]] const std::string& superPopulation() const { return super_population_; }
  [[nodiscard]] double minAlleleFrequency() const { return allele_frequency_min_; }
  [[nodiscard]] double maxAlleleFrequency() const { return allele_frequency_max_; }
  [[nodiscard]] std::shared_ptr<const ContigVariant> unphasedContigPtr() const { return unphased_contig_ptr_; }
  [[nodiscard]] VariantDatabaseSource variantSource() const { return variant_source_; }

  void lociiType(LociiType type) { locii_type_ = type; }
  void lowerOffset(ContigOffset_t lower) { lower_offset_ = lower; }
  void upperOffset(ContigOffset_t upper) { upper_offset_ = upper; }
  void lociiSpacing(size_t spacing) { spacing_ = spacing; };
  void lociiCount (size_t count) { locii_count = count; }
  void superPopulation(const std::string& super_population){ super_population_ = super_population; }
  void minAlleleFrequency(double min_AF) { allele_frequency_min_ = std::clamp(min_AF, 0.0, 1.0); }
  void maxAlleleFrequency(double max_AF) { allele_frequency_max_ = std::clamp(max_AF, 0.0, 1.0); }
  void variantSource(VariantDatabaseSource variant_source) { variant_source_ = variant_source; }

private:

  std::shared_ptr<const ContigVariant> unphased_contig_ptr_;
  LociiType locii_type_{LociiType::LOCII_FROM_COUNT};
  ContigOffset_t lower_offset_{0};
  ContigOffset_t upper_offset_{1000000000};
  size_t spacing_{1000};
  size_t locii_count{1000};
  std::string super_population_{"EUR"};
  double allele_frequency_min_{0.0};
  double allele_frequency_max_{1.0};
  VariantDatabaseSource variant_source_{VariantDatabaseSource::GNOMAD2_1};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object uses a range of selection criteria to generate a vector of locii offsets used to select locii for
// inbreeding sampling.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class RetrieveLociiVector {

public:

  RetrieveLociiVector() =delete;
  ~RetrieveLociiVector() = delete;

  static std::vector<ContigOffset_t> getLociiVector(const LociiVectorArguments& arguments);

private:


  static std::vector<ContigOffset_t> getLociiFromTo(std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                             ContigOffset_t lower_offset,
                                             ContigOffset_t upper_offset,
                                             size_t spacing,
                                             const std::string& super_population,
                                             double allele_frequency_min,
                                             double allele_frequency_max,
                                             VariantDatabaseSource variant_source);

  static std::vector<ContigOffset_t> getLociiCount( std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                             ContigOffset_t start_offset,
                                             size_t locii_count,
                                             size_t spacing,
                                             const std::string& super_population,
                                             double allele_frequency_min,
                                             double allele_frequency_max,
                                             VariantDatabaseSource variant_source);

  static std::vector<AlleleFreqVector> getAllelesFromTo( std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                                         ContigOffset_t lower_offset,
                                                         ContigOffset_t upper_offset,
                                                         size_t spacing,
                                                         const std::string& super_population,
                                                         double allele_frequency_min,
                                                         double allele_frequency_max,
                                                         VariantDatabaseSource variant_source);

  static std::vector<AlleleFreqVector> getAllelesCount( std::shared_ptr<const ContigVariant> unphased_contig_ptr,
                                                         ContigOffset_t lower_offset,
                                                         size_t locii_count,
                                                         size_t spacing,
                                                         const std::string& super_population,
                                                         double allele_frequency_min,
                                                         double allele_frequency_max,
                                                         VariantDatabaseSource variant_source);

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object actually uses the locii vector generated above to package the selected locii into a ContigLocusMap
// for presentation to the inbreeding analytics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Allele locus lists indexed by superpopulation
using LocusMap = std::map<std::string, std::shared_ptr<const ContigVariant>>;
// LocusMaps indexed by contig id.
using ContigLocusMap = std::map<ContigId_t, LocusMap>;

// Just a namespace.
class InbreedSampling  {

public:

  InbreedSampling() = delete;
  ~InbreedSampling() = delete;


  // Uses the defined contigs in the unphased population to create a contig map of population locii.
  [[nodiscard]] static ContigLocusMap getPopulationLocusMap(  std::shared_ptr<const UnphasedPopulation> population_ptr,
                                                              double min_af,
                                                              double max_af,
                                                              ContigOffset_t locii_spacing,
                                                              ContigOffset_t upper_offset = DEFAULT_UPPER,
                                                              ContigOffset_t lower_offset = DEFAULT_LOWER);

  constexpr static const ContigOffset_t DEFAULT_UPPER = 1000000000;
  constexpr static const ContigOffset_t DEFAULT_LOWER = 0;

private:


  // Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
  // and at a specified frequency for the super population. Used as a template for calculating
  // the inbreeding coefficient and sample relatedness
  using LocusReturnPair = std::pair<std::string, std::shared_ptr<const ContigVariant>>;
  [[nodiscard]] static LocusReturnPair getLocusList( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                     const ContigId_t& contig_id,
                                                     const std::string& super_population,
                                                     double min_af,
                                                     double max_af,
                                                     ContigOffset_t locii_spacing,
                                                     ContigOffset_t upper_offset,
                                                     ContigOffset_t lower_offset);

  // Generate a list of locii to sample a population for the inbreeding coefficient.
  [[nodiscard]] static LocusMap getPopulationLocus(std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const ContigId_t& contig_id,
                                                   double min_af,
                                                   double max_af,
                                                   ContigOffset_t locii_spacing,
                                                   ContigOffset_t upper_offset,
                                                   ContigOffset_t lower_offset);



};





} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_INBREED_AUX_H
