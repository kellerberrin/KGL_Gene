//
// Created by kellerberrin on 22/11/20.
//

#ifndef KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H
#define KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H

#include "kgl_variant_db_phased.h"
#include "kgl_analysis_mutation_inbreed_freqdb.h"

#include <memory>
#include <map>

namespace kellerberrin::genome {   //  organization::project level namespace


// Just a namespace.
class InbreedSynthetic  {

public:

  InbreedSynthetic() = delete;
  ~InbreedSynthetic() = delete;


  // Create a synthetic population with known inbreeding characteristics
  // Used to test and calibrate the developed inbreeding algorithms.
  [[nodiscard]] static std::shared_ptr<const DiploidPopulation>
  generateSyntheticPopulation( double lower_inbreeding,
                               double upper_inbreeding,
                               double step_inbreeding,
                               const std::string& super_population,
                               const ContigVariant& locus_list);


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



#endif //KGL_KGL_ANALYSIS_MUTATION_SYNTHETIC_H
