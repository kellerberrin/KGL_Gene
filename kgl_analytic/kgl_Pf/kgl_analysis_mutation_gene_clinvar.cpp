//
// Created by kellerberrin on 31/1/21.
//

#include "kgl_analysis_mutation_gene_clinvar.h"
#include "kgl_variant_filter.h"


namespace kgl = kellerberrin::genome;




bool kgl::VariantPhaseStats::phaseStatistics(const std::shared_ptr<const ContigDB>& subject_variants) {


  auto male_variants = subject_variants->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_B));
  size_t male_phase_variants = male_variants->variantCount();
  if (male_phase_variants > 0) {

    ++phase_male_;

  }

  auto female_variants = subject_variants->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_A));
  size_t female_phase_variants = female_variants->variantCount();
  if (female_phase_variants > 0) {

    ++phase_female_;

  }

  if (male_phase_variants > 0 and female_phase_variants > 0) {

    ++phase_hom_;

  }

  if (male_phase_variants > 0 or female_phase_variants > 0) {

    ++phase_either_;

  }

  return true;

}


void kgl::VariantPhaseStats::writeHeader(std::ostream& out_file, char output_delimiter) {

  out_file << "PhCount" << output_delimiter
           << "PhMale" << output_delimiter
           << "PhFemale" << output_delimiter
           << "PhHom";

}
