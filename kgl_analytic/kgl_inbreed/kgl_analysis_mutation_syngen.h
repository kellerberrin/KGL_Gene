//
// Created by kellerberrin on 30/11/20.
//

#ifndef KGL_ANALYSIS_MUTATION_SYNGEN_H
#define KGL_ANALYSIS_MUTATION_SYNGEN_H


#include "kgl_analysis_mutation_inbreed_freqdb.h"
#include "kgl_analysis_mutation_inbreed_locus.h"


#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Generate synthetic genomes with known inbreeding coefficients. These are generated stochastically using
// minor allele frequencies to imitate real genomes as closely as possible. Used to test inbreeding analytics.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Just a namespace.
class InbreedSynthetic  {

public:

  InbreedSynthetic() = delete;
  ~InbreedSynthetic() = delete;


  // Create a synthetic population with known inbreeding characteristics
  // Used to test and calibrate the developed inbreeding algorithms.
  [[nodiscard]] static std::shared_ptr<const PopulationDB>
  generateSyntheticPopulation( double lower_inbreeding,
                               double upper_inbreeding,
                               double step_inbreeding,
                               const std::string& super_population,
                               const ContigDB& locus_list,
                               const LociiVectorArguments& arguments);


  // Recover a synthetic inbreeding coefficient from the synthetic genome id.
  [[nodiscard]] static std::pair<bool, double> generateInbreeding(const GenomeId_t& genome_id);


  // Synthetic genome constant
  constexpr static const double SYNTHETIC_GENOME = 1000000; // Used to create the synthetic genome id.

private:

  // Generate an inbreeding encoded synthetic genome
  [[nodiscard]] static GenomeId_t generateSyntheticGenomeId( double inbreeding,
                                                             const std::string& super_population,
                                                             size_t counter);

};



} // namespace






#endif //KGL_KGL_ANALYSIS_MUTATION_SYNGEN_H
