//
// Created by kellerberrin on 24/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_AUX_H
#define KGL_ANALYSIS_MUTATION_INBREED_AUX_H

#include "kgl_analysis_inbreed_freq.h"

#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace

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
// This object uses a range of selection criteria to generate a vector of locii offsets used to select locii for
// inbreeding sampling.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class RetrieveLociiVector {

public:

  RetrieveLociiVector() =delete;
  ~RetrieveLociiVector() = delete;

  static std::vector<ContigOffset_t> getLociiFromTo(std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                    const std::string& super_population,
                                                    const LociiVectorArguments& arguments);

  static std::vector<ContigOffset_t> getLociiCount(std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                   const std::string& super_population,
                                                   const LociiVectorArguments& arguments);

private:


  static std::vector<AlleleFreqVector> getAllelesFromTo( std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                         const std::string& super_population,
                                                         const LociiVectorArguments& arguments);

  static std::vector<AlleleFreqVector> getAllelesCount( std::shared_ptr<const ContigDB> unphased_contig_ptr,
                                                        const std::string& super_population,
                                                        const LociiVectorArguments& arguments);

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object actually uses the locii vector generated above to package the selected locii into a ContigLocusMap
// for presentation to the inbreeding analytics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Allele locus lists indexed by superpopulation
using LocusMap = std::map<std::string, std::shared_ptr<const ContigDB>>;
// LocusMaps indexed by contig id.
using ContigLocusMap = std::map<ContigId_t, LocusMap>;

// Just a namespace.
class InbreedSampling  {

public:

  InbreedSampling() = delete;
  ~InbreedSampling() = delete;


  // Uses the defined contigs in the unphased population to create a contig map of population locii.
  [[nodiscard]] static ContigLocusMap getPopulationLocusMap(  std::shared_ptr<const PopulationDB> population_ptr,
                                                              const LociiVectorArguments& locii_args);

private:


  // Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
  // and at a specified frequency for the super population. Used as a template for calculating
  // the inbreeding coefficient and sample relatedness
  using LocusReturnPair = std::pair<std::string, std::shared_ptr<const ContigDB>>;
  [[nodiscard]] static LocusReturnPair getLocusList( std::shared_ptr<const PopulationDB> unphased_ptr,
                                                     const ContigId_t& contig_id,
                                                     const std::string& super_population,
                                                     const LociiVectorArguments& locii_args);

  // Generate a list of locii to sample a population for the inbreeding coefficient.
  [[nodiscard]] static LocusMap getPopulationLocus(std::shared_ptr<const PopulationDB> unphased_ptr,
                                                   const ContigId_t& contig_id,
                                                   const LociiVectorArguments& locii_args);

};





} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_INBREED_AUX_H
