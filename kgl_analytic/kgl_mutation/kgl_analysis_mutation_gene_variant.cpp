//
// Created by kellerberrin on 4/3/21.
//

#include "kgl_analysis_mutation_gene_variant.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"


namespace kgl = kellerberrin::genome;


void kgl::GeneVariants::ProcessVariantStats(const GenomeId_t& genome_id,
                                            const std::shared_ptr<const ContigDB>& span_variant_ptr,
                                            const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                            const std::shared_ptr<const GenomePEDData>& ped_data) {

  std::map<std::string, size_t> variant_distribution;
  // Variant statistics.
  ++genome_count_;

  // Get phased VEP info.
  VepInfo lof_het_hom = geneSpanVep( span_variant_ptr, unphased_population_ptr);

  if (lof_het_hom.all_lof > 0) {

    ++all_lof_;
    updateEthnicity().pedAnalysis(genome_id, 1, ped_data);

  }


  if (lof_het_hom.female_phase_lof > 0) {

    ++female_lof_;

  }

  if (lof_het_hom.male_phase_lof > 0) {

    ++male_lof_;

  }

  // Loss of Function in both chromosomes, Mendelian genetics.
  if (lof_het_hom.female_phase_lof > 0 and lof_het_hom.male_phase_lof > 0) {

    ++hom_lof_;

  }


  if (lof_het_hom.female_high_effect > 0) {

    ++female_high_effect_;

  }

  if (lof_het_hom.male_high_effect > 0) {

    ++male_high_effect_;

  }

  // High Impact in both chromosomes, Mendelian genetics.
  if (lof_het_hom.male_high_effect > 0 and lof_het_hom.female_high_effect > 0) {

    ++hom_high_effect_;

  }

  size_t variant_count = span_variant_ptr->variantCount();

  if (variant_count > 0) {

    ++genome_variant_;
    span_variant_count_ += variant_count;

  }


  double indel{0.0}, transition{0.0}, transversion{0.0};

  for (auto const& [offset, offset_ptr] : span_variant_ptr->getMap()) {

    const OffsetDBArray& variant_array = offset_ptr->getVariantArray();

    for (auto const& variant_ptr : variant_array) {


      if (variant_ptr->phaseId() == VariantPhase::DIPLOID_PHASE_A) {

        ++female_phase_;

      } else {

        ++male_phase_;

      }

      switch (variant_ptr->variantType()) {

        case VariantType::INDEL:
          indel += 1.0;
          break;

        case VariantType::TRANSITION:
          transition += 1.0;
          break;

        case VariantType::TRANSVERSION:
          transversion += 1.0;
          break;

      }

      auto find_result = variant_distribution.find(variant_ptr->variantPhaseHash());
      if (find_result == variant_distribution.end()) {

        auto[it, insert_result] = variant_distribution.try_emplace(variant_ptr->variantPhaseHash(), 1);
        if (not insert_result) {

          ExecEnv::log().error("GenomeMutation::variantAnalysis; cannot insert (duplicate) variant with hash: {}",
                               variant_ptr->variantPhaseHash());

        }

      } else {

        auto&[variant_ident, count] = *find_result;
        ++count;

      }

    } //for variant

  } //for offset

  double sum = indel + transition + transversion;
  if (sum >= 1.0) {

    indel = indel / sum;
    transition = transition / sum;
    transversion = transversion / sum;

  } else {

    indel = 0.0;
    transition = 0.0;
    transversion = 0.0;

  }

  unique_variants_ += variant_distribution.size();
  ++variant_count_;

  for (auto const& [offset, offset_ptr] : span_variant_ptr->getMap()) {

    if (offset_ptr->getVariantArray().size() == 1) {

      ++heterozygous_;

    } else if (offset_ptr->getVariantArray().size() == 2) {

      auto const& offset_array = offset_ptr->getVariantArray();
      if (offset_array.front()->variantHash() == offset_array.back()->variantHash()) {

        ++homozygous_;

      }

    }

  }

}



kgl::VepInfo kgl::GeneVariants::geneSpanVep( const std::shared_ptr<const ContigDB>& span_contig,
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

  auto phase_A_variants = span_contig->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_A));
  auto found_unphased_A = unphased_contig->findContig(phase_A_variants);


  auto phase_B_variants = span_contig->filterVariants(PhaseFilter(VariantPhase::DIPLOID_PHASE_B));
  auto found_unphased_B = unphased_contig->findContig(phase_B_variants);

  auto found_all_unphased = unphased_contig->findContig(span_contig);

  vep_info.all_lof = VepCount(found_all_unphased, LOF_VEP_FIELD_, LOF_HC_VALUE_);
  vep_info.female_phase_lof = VepCount(found_unphased_A, LOF_VEP_FIELD_, LOF_HC_VALUE_);
  vep_info.male_phase_lof = VepCount(found_unphased_B, LOF_VEP_FIELD_, LOF_HC_VALUE_);

  VepSubFieldValues vep_field(IMPACT_VEP_FIELD_);

  vep_field.getContigValues(found_unphased_A);

  auto result = vep_field.getMap().find(IMPACT_HIGH_VALUE_);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.female_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE_);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.female_moderate_effect = count;

  }

  vep_field.getContigValues(found_unphased_B);

  result = vep_field.getMap().find(IMPACT_HIGH_VALUE_);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.male_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE_);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.male_moderate_effect = count;

  }

  return vep_info;

}


size_t kgl::GeneVariants::VepCount( const std::shared_ptr<const ContigDB>& vep_contig,
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

