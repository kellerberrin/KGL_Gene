//
// Created by kellerberrin on 27/1/21.
//


#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_analysis_mutation_gene.h"
#include "kgl_variant_mutation.h"

#include <fstream>


namespace kgl = kellerberrin::genome;




std::shared_ptr<const kgl::ContigDB> kgl::GenomeMutation::getGeneSpan(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                      const GeneCharacteristic& gene_char) {

  return contig_ptr->subset(gene_char.gene_begin, gene_char.gene_end);

}




kgl::VepInfo kgl::GenomeMutation::geneSpanVep( const std::shared_ptr<const ContigDB>& span_contig,
                                               const std::shared_ptr<const PopulationDB>& unphased_population_ptr) {

  VepInfo vep_info;

  if (unphased_population_ptr->getMap().size() != 1) {

    ExecEnv::log().error("GenomeMutation::geneSpanVep; expected unphased population to have 1 genome, size if: {}", unphased_population_ptr->getMap().size());
    return vep_info;

  }

  auto [genomne_id, genome_ptr] = *(unphased_population_ptr->getMap().begin());

  auto contig_opt = genome_ptr->getContig(span_contig->contigId());

  if (not contig_opt) {

    return vep_info;

  }

  auto unphased_contig = contig_opt.value();

  auto phase_A_variants = span_contig->filterVariants(PhaseFilter(VariantSequence::DIPLOID_PHASE_A));
  auto found_unphased_A = unphased_contig->findContig(phase_A_variants);


  auto phase_B_variants = span_contig->filterVariants(PhaseFilter(VariantSequence::DIPLOID_PHASE_B));
  auto found_unphased_B = unphased_contig->findContig(phase_B_variants);

  vep_info.female_lof = VepCount(found_unphased_A, LOF_VEP_FIELD, LOF_HC_VALUE);
  vep_info.male_lof = VepCount(found_unphased_B, LOF_VEP_FIELD, LOF_HC_VALUE);

  VepSubFieldValues vep_field(IMPACT_VEP_FIELD);

  vep_field.getContigValues(phase_A_variants);

  if (vep_field.getMap().size() > 0) ExecEnv::log().info("GenomeMutation::geneSpanVep; IMPACT entries: {}", vep_field.getMap().size());

  auto result = vep_field.getMap().find(IMPACT_HIGH_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.female_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.female_moderate_effect = count;

  }

  vep_field.getContigValues(phase_B_variants);

  result = vep_field.getMap().find(IMPACT_HIGH_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.male_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.male_moderate_effect = count;

  }

  return vep_info;

}


size_t kgl::GenomeMutation::VepCount( const std::shared_ptr<const ContigDB>& vep_contig,
                                      const std::string& vep_field_ident,
                                      const std::string& vep_field_value) {

  VepSubFieldValues vep_field(vep_field_ident);
  vep_field.getContigValues(vep_contig);

  auto result = vep_field.getMap().find(vep_field_value);
  if (result != vep_field.getMap().end()) {

    auto [field_value, count] = *result;
    return count;

  }

  return 0;

}



void kgl::GenomeMutation::updatePopulations(const std::shared_ptr<const GenomePEDData>& ped_data) {

  // Generate the template populations.
  // Insert the templates into the gene records.
  for (auto& gene_record : gene_vector_) {

    gene_record.clinvar.results.updatePopulations(ped_data);

  }

}
