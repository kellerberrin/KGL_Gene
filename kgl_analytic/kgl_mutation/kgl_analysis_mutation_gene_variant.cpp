//
// Created by kellerberrin on 4/3/21.
//

#include "kgl_analysis_mutation_gene_variant.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kel_distribution.h"


namespace kgl = kellerberrin::genome;


void kgl::GeneVariants::processVariantStats(const GenomeId_t& genome_id,
                                            const std::shared_ptr<const ContigDB>& span_variant_ptr,
                                            const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                            const std::shared_ptr<const GenomePEDData>& ped_data) {

  std::map<std::string, size_t> variant_distribution;
  // Variant statistics.
  ++genome_count_;

  // Get phased VEP info.
  VepInfo vep_info = geneSpanVep(span_variant_ptr, unphased_population_ptr);

  if (vep_info.all_lof > 0) {

    ++all_lof_;
    updateLofEthnicity().pedAnalysis(genome_id, 1, ped_data);

  }


  if (vep_info.female_phase_lof > 0) {

    ++female_lof_;

  }

  if (vep_info.male_phase_lof > 0) {

    ++male_lof_;

  }

  // Loss of Function in both chromosomes, Mendelian genetics.
  if (vep_info.hom_lof > 0) {

    ++hom_lof_;

  }

  if (vep_info.all_high_effect > 0) {

    ++all_high_effect_;
    updateHighEthnicity().pedAnalysis(genome_id, 1, ped_data);

  }

  // High Impact in both chromosomes, Mendelian genetics.
  if (vep_info.hom_high_effect > 0) {

    ++hom_high_effect_;

  }

  if (vep_info.all_moderate_effect > 0) {

    ++all_moderate_effect_;
    updateModerateEthnicity().pedAnalysis(genome_id, 1, ped_data);

  }

  // High Impact in both chromosomes, Mendelian genetics.
  if (vep_info.hom_moderate_effect > 0) {

    ++hom_moderate_effect_;

  }

  size_t variant_count = span_variant_ptr->variantCount();

  if (variant_count > 0) {

    ++genome_variant_;
    span_variant_count_ += variant_count;

  }

  for (auto const& [offset, offset_ptr] : span_variant_ptr->getMap()) {

    const OffsetDBArray& variant_array = offset_ptr->getVariantArray();

    for (auto const& variant_ptr : variant_array) {


      if (variant_ptr->phaseId() == VariantPhase::DIPLOID_PHASE_A) {

        ++female_phase_;

      } else {

        ++male_phase_;

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

  std::shared_ptr<const ContigDB> phase_hom_variants = found_unphased_B->setIntersection(*found_unphased_A, VariantEquality::UNPHASED);

  vep_info.all_lof = vepCount(found_all_unphased, LOF_VEP_FIELD_, LOF_HC_VALUE_);
  vep_info.female_phase_lof = vepCount(found_unphased_A, LOF_VEP_FIELD_, LOF_HC_VALUE_);
  vep_info.male_phase_lof = vepCount(found_unphased_B, LOF_VEP_FIELD_, LOF_HC_VALUE_);
  vep_info.hom_lof = vepCount(phase_hom_variants, LOF_VEP_FIELD_, LOF_HC_VALUE_);


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

  vep_field.getContigValues(phase_hom_variants);

  result = vep_field.getMap().find(IMPACT_HIGH_VALUE_);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.hom_high_effect = count;

  }

  result = vep_field.getMap().find(IMPACT_MODERATE_VALUE_);
  if (result != vep_field.getMap().end()) {

    auto [field_ident, count] = *result;
    vep_info.hom_moderate_effect = count;

  }

  vep_info.all_high_effect = vep_info.male_high_effect + vep_info.female_high_effect;
  vep_info.all_moderate_effect = vep_info.male_moderate_effect + vep_info.female_moderate_effect;


  return vep_info;

}


size_t kgl::GeneVariants::vepCount( const std::shared_ptr<const ContigDB>& vep_contig,
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


bool kgl::GeneVariants::processSummaryStatistics( const std::shared_ptr<const PopulationDB> &population_ptr,
                                                  const GeneEthnicitySex& ethnic_statistics,
                                                  const std::string& gene) {


  bool tail_result{true};
  size_t total_success = ethnic_lof_.total() + ethnic_high_.total() + ethnic_moderate_.total();
  size_t population_size = population_ptr->getMap().size() * 3;

  if (not ethnic_lof_.auditTotals()) {

    ExecEnv::log().error("GeneVariants::processSummaryStatistics; problem with Vep LoF ethnic statistics totals");
    tail_result = false;

  }

  if (not ethnic_high_.auditTotals()) {

    ExecEnv::log().error("GeneVariants::processSummaryStatistics; problem with Vep High ethnic statistics totals");
    tail_result = false;

  }

  if (not ethnic_moderate_.auditTotals()) {

    ExecEnv::log().error("GeneVariants::processSummaryStatistics; problem with Vep Moderate ethnic statistics totals");
    tail_result = false;

  }

  for (auto const& [super_population, sample_size] : ethnic_statistics.superPopulation()) {

    double upper_tail{0.0};
    double lower_tail{0.0};

    if (total_success > 0) {

      size_t lof_count = ethnic_lof_.superPopulationCount(super_population);
      size_t high_count = ethnic_high_.superPopulationCount(super_population);
      size_t mod_count = ethnic_moderate_.superPopulationCount(super_population);

      size_t pop_success = lof_count + high_count + mod_count;

      size_t variant_sample_size = sample_size * 3;

      HypergeometricDistribution hyper_tail_stats(total_success, variant_sample_size, population_size);

      if (pop_success > hyper_tail_stats.upperSuccesses_k()) {

        ExecEnv::log().warn( "GeneVariants::processSummaryStatistics; Gene: {}, hypergeometric parameters, K: {}, N: {}, k: {}, n: {}",
                             gene, total_success, population_size, pop_success, sample_size);
        ExecEnv::log().warn( "GeneVariants::processSummaryStatistics; Gene: {}, hypergeometric parameters, lof_count: {}, high_count: {}, moderate_count: {}",
                             gene, lof_count, high_count, mod_count);
        tail_result = false;

      }

      upper_tail = hyper_tail_stats.upperSingleTailTest(pop_success);
      lower_tail = hyper_tail_stats.lowerSingleTailTest(pop_success);


    } // success > 0

    auto upper_find = upper_tail_.find(super_population);
    if (upper_find != upper_tail_.end()) { // Found

      auto& [population, tail_value] = *upper_find;
      if (tail_value != 0.0 and total_success > 0) {

        ExecEnv::log().error("GeneVariants::processSummaryStatistics; found unexpected non-zero upper tail statistic: {}, population: {}, ",
                             tail_value,  population);
        tail_value = upper_tail;
        tail_result = false;

      } else if (tail_value == 0.0 and total_success != 0) {

        tail_value = lower_tail;

      }

    } else { // Not found.

      upper_tail_.emplace(super_population, upper_tail);

    }

    auto lower_find = lower_tail_.find(super_population);
    if (lower_find != lower_tail_.end()) {

      auto& [population, tail_value] = *lower_find;
      if (tail_value != 0.0 and total_success != 0) {

        ExecEnv::log().error("GeneVariants::processSummaryStatistics; found unexpected non-zero lower tail statistic: {}, population: {}",
                             tail_value, population);
        tail_value = lower_tail;
        tail_result = false;

      } else if (tail_value == 0.0 and total_success != 0) {

        tail_value = lower_tail;

      }

    } else {

      lower_tail_.emplace(super_population, lower_tail);

    }

  }

  return tail_result;

}

void kgl::GeneVariants::initializeSummaryStatistics( const GeneEthnicitySex& ethnic_statistics) {

  if (upper_tail_.size() == ethnic_statistics.superPopulation().size()
      and lower_tail_.size() == ethnic_statistics.superPopulation().size()) {

    return;

  }

  upper_tail_.clear();
  lower_tail_.clear();

  for (auto const& [super_population, sample_size] : ethnic_statistics.superPopulation()) {

    upper_tail_.emplace(super_population, 0.0);
    lower_tail_.emplace(super_population, 0.0);

  }

}
