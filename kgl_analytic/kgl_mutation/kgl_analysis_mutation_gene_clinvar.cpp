//
// Created by kellerberrin on 31/1/21.
//

#include "kgl_analysis_mutation_gene_clinvar.h"
#include "kgl_analysis_mutation_clinvar.h"
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::GeneClinvar::writeOutput(  const std::shared_ptr<const GenomePEDData>& ped_data,
                                         std::ostream& out_file,
                                         char output_delimiter) const {

  out_file << getPhase().phaseEither() << output_delimiter
           << getPhase().phaseMale() << output_delimiter
           << getPhase().phaseFemale() << output_delimiter
           << getPhase().phaseHomozygous() << output_delimiter;

  std::string concat_desc{"\""};
  for (auto const& desc : getClinvarDesc()) {

    if (desc != *getClinvarDesc().begin()) {

      concat_desc += CONCAT_TOKEN_;

    }

    concat_desc += desc;

  }
  concat_desc += "\"";
  out_file << concat_desc << output_delimiter;

  getEthnicity().writeOutput(ped_data, out_file, output_delimiter);

}


void kgl::GeneClinvar::writeHeader( const std::shared_ptr<const GenomePEDData>& ped_data,
                                    std::ostream& out_file,
                                    char output_delimiter) const {

  VariantPhaseStats::writeHeader(out_file, output_delimiter);

  out_file << output_delimiter
           << "CLV_Desc" << output_delimiter;

  getEthnicity().writeHeader(ped_data, out_file, output_delimiter);

}



void kgl::GeneClinvar::processClinvar( const GenomeId_t& genome_id,
                                       const std::shared_ptr<const ContigDB>& subject_variants,
                                       const std::shared_ptr<const ContigDB>& clinvar_contig,
                                       const std::shared_ptr<const GenomePEDData>& ped_data) {


  auto subject_clinvar = AnalyzeClinvar::findClinvar(subject_variants, clinvar_contig);
  auto info_vector = AnalyzeClinvar::clinvarInfo(subject_clinvar);
  size_t clinvar_variants = subject_clinvar->variantCount();
  if (clinvar_variants > 0) {

    ++updatePhase().updatePhaseEither();

  }

  for (auto const& record : info_vector) {

    // Save any unique descriptions.
    updateClinvarDesc().insert(record.clndn);

  }

  auto male_variants = subject_variants->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_B));
  auto male_clinvar = AnalyzeClinvar::findClinvar(male_variants, clinvar_contig);
  size_t male_phase_variants = male_clinvar->variantCount();
  if (male_phase_variants > 0) {

    ++updatePhase().updatePhaseMale();

  }

  auto female_variants = subject_variants->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_A));
  auto female_clinvar = AnalyzeClinvar::findClinvar(female_variants, clinvar_contig);
  size_t female_phase_variants = female_clinvar->variantCount();
  if (female_phase_variants > 0) {

    ++updatePhase().updatePhaseFemale();

  }

  if (male_phase_variants > 0 and female_phase_variants > 0) {

    ++updatePhase().updatePhaseHomozygous();

  }

  size_t count = clinvar_variants > 0 ? 1 : 0;
  updateEthnicity().pedAnalysis(genome_id, count, ped_data);

}


