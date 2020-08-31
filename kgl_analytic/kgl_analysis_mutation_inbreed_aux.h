//
// Created by kellerberrin on 24/8/20.
//

#ifndef KGL_ANALYSIS_MUTATION_INBREED_AUX_H
#define KGL_ANALYSIS_MUTATION_INBREED_AUX_H

#include "kgl_variant_db_phased.h"

#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace


// Allele locus lists indexed by superpopulation
using LocusMap = std::map<std::string, std::shared_ptr<const ContigVariant>>;
// LocusMaps indexed by contig id.
using ContigLocusMap = std::map<ContigId_t, LocusMap>;

// Just a namespace.
class InbreedSampling  {

public:

  InbreedSampling() = delete;
  ~InbreedSampling() = delete;

  // Generate a list of locii to sample a population for the inbreeding coefficient.
  [[nodiscard]] static LocusMap getPopulationLocus(std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                   const ContigId_t& contig_id,
                                                   double min_af = MIN_MAF_,
                                                   double max_af = MAX_MAF_,
                                                   ContigOffset_t locii_spacing = VARIANT_SPACING_);

  // Uses the defined contigs in the unphased population to create a contig map of population locuii.
  [[nodiscard]] static ContigLocusMap getPopulationLocusMap(  std::shared_ptr<const UnphasedPopulation> population_ptr,
                                                              double min_af,
                                                              double max_af,
                                                              ContigOffset_t locii_spacing);
  // Use a super population code to lookup a corresponding AF field.
  [[nodiscard]] static std::string lookupSuperPopulationField(const std::string& super_population);

  // Inverse of the above.
  [[nodiscard]] static std::string inverseSuperPopulationField(const std::string& super_population);

  // Get a scalar floating Info field.
  [[nodiscard]] static std::pair<bool, double> processFloatField( const Variant& variant,
                                                                  const std::string& field_name);

  // Create a synthetic population with known inbreeding characteristics
  // Used to test and calibrate the developed inbreeding algorithms.
  [[nodiscard]] static std::shared_ptr<const DiploidPopulation>
    generateSyntheticPopulation( double lower_inbreeding,
                                 double upper_inbreeding,
                                 double step_inbreeding,
                                 const std::string& super_population,
                                 const ContigVariant& locus_list);

  constexpr static const double SELECTION_PERCENT_ = 0.20; // Select 20 percent of the locus variants.

  // Recover a synthetic inbreeding coefficient from the synthetic genome id.
  [[nodiscard]] static std::pair<bool, double> generateInbreeding(const GenomeId_t& genome_id);

  // Generate an inbreeding encoded synthetic genome
  [[nodiscard]] static GenomeId_t generateSyntheticGenomeId( double inbreeding,
                                                             const std::string& super_population,
                                                             size_t counter);


  constexpr static const double MIN_MAF_ = 0.10; // Set the minimum MAF
  constexpr static const double MAX_MAF_ = 1.05; // Set the maximum MAF
  constexpr static const ContigOffset_t VARIANT_SPACING_ = 1000;  // Set the variant spacing.
  // Synthetic genome constant
  constexpr static const double SYNTHETIC_GENOME = 1000000; // Used to create the synthetic genome id.

private:


  // Get a list of potential allele locus with a specified spacing to minimise linkage dis-equilibrium
  // and at a specified frequency for the super population. Used as a template for calculating
  // the inbreeding coefficient and sample relatedness
  using SuperPopPair = std::pair<const char*, const char*>;
  using LocusReturnPair = std::pair<std::string, std::shared_ptr<const ContigVariant>>;
  [[nodiscard]] static LocusReturnPair getLocusList( std::shared_ptr<const UnphasedPopulation> unphased_ptr,
                                                     const ContigId_t& contig_id,
                                                     ContigOffset_t spacing,
                                                     const SuperPopPair& super_pop_pair,
                                                     double min_frequency,
                                                     double max_frequency);

  // The Info field identifiers for allele frequency for Gnomad 2.1 and 3.0
  // First is the super populations defined in the 1000 Genomes.
  // Second is the population lookup in the Gnomad Info data.
  constexpr static const SuperPopPair SUPER_POP_AFR_GNOMAD_ {"AFR", "AF_afr"} ;  // African
  constexpr static const SuperPopPair SUPER_POP_AMR_GNOMAD_ {"AMR", "AF_amr"};  // American
  constexpr static const SuperPopPair SUPER_POP_EAS_GNOMAD_ {"EAS", "AF_eas"};  // East Asian
  constexpr static const SuperPopPair SUPER_POP_EUR_GNOMAD_ { "EUR", "AF_nfe"};  // European
  constexpr static const SuperPopPair SUPER_POP_SAS_GNOMAD_ { "SAS", "AF"};  // South Asian
  // Use a super population code to lookup a corresponding AF field.

  inline static const std::vector<SuperPopPair> super_pop_vec { SUPER_POP_AFR_GNOMAD_,
                                                                SUPER_POP_AMR_GNOMAD_,
                                                                SUPER_POP_EAS_GNOMAD_,
                                                                SUPER_POP_EUR_GNOMAD_,
                                                                SUPER_POP_SAS_GNOMAD_};


};





} // namespace

#endif //KGL_KGL_ANALYSIS_MUTATION_INBREED_AUX_H
