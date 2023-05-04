//
// Created by kellerberrin on 4/3/21.
//

#include "kgl_analysis_mutation_gene_variant.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_factory_vcf_evidence_vep.h"
#include "kgl_variant_filter_db.h"
#include "kel_distribution.h"


namespace kgl = kellerberrin::genome;

void kgl::GeneVariants::initializeEthnic(const std::shared_ptr<const HsGenomeAux>& genome_aux_data) {

  ethnic_lof_.setDisplay("LOF_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));
  ethnic_high_.setDisplay("HIGH_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));
  ethnic_moderate_.setDisplay("MOD_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));
  ethnic_citation_.setDisplay("CITE_", (GeneEthnicitySex::DISPLAY_SEX_FLAG | GeneEthnicitySex::DISPLAY_SUPER_POP_FLAG));

  ethnic_lof_.updatePopulations(genome_aux_data);
  ethnic_high_.updatePopulations(genome_aux_data);
  ethnic_moderate_.updatePopulations(genome_aux_data);
  ethnic_citation_.updatePopulations(genome_aux_data);

}


void kgl::GeneVariants::processVariantStats(const GenomeId_t& genome_id,
                                            const std::shared_ptr<const ContigDB>& span_variant_ptr,
                                            const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                            const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                            const std::shared_ptr<const CitationResource>& allele_citation_ptr) {

  // Variant statistics.
  ++genome_count_;

  // Get phased VEP info.
  VepInfo vep_info = geneSpanVep(span_variant_ptr, unphased_population_ptr);

  if (vep_info.all_lof > 0) {

    ++all_lof_;
    ethnic_lof_.genomeAnalysis(genome_id, 1, genome_aux_data);

  }

  // Loss of Function in both chromosomes, Mendelian genetics.
  if (vep_info.hom_lof > 0) {

    ++hom_lof_;

  }

  if (vep_info.all_high_effect > 0) {

    ++all_high_effect_;
    ethnic_high_.genomeAnalysis(genome_id, 1, genome_aux_data);

  }

  // High Impact in both chromosomes, Mendelian genetics.
  if (vep_info.hom_high_effect > 0) {

    ++hom_high_effect_;

  }

  if (vep_info.all_moderate_effect > 0) {

    ++all_moderate_effect_;
    ethnic_moderate_.genomeAnalysis(genome_id, 1, genome_aux_data);

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

      unique_variants_.insert(variant_ptr->HGVS());
// count variants with citations.
      if (not variant_ptr->identifier().empty()) {

        auto find_result = allele_citation_ptr->alleleIndexedCitations().find(variant_ptr->identifier());
        if (find_result != allele_citation_ptr->alleleIndexedCitations().end()) {

          auto const& [rsid, citations] = *find_result;
          for (auto const& cite : citations) {

            unique_citations_.insert(cite);

          }
          // Increment if any citation
          ++citation_count_;
          // Ethnic profile
          ethnic_citation_.genomeAnalysis(genome_id, 1, genome_aux_data);

        }

      }

    } //for variant

  } //for offset

  ++variant_count_;

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

  auto found_all_unphased = unphased_contig->findContig(span_contig);
  std::shared_ptr<const ContigDB> found_hom_variants = span_contig->viewFilter(HomozygousFilter());

  vep_info.all_lof = vepCount(found_all_unphased, LOF_VEP_FIELD_, LOF_HC_VALUE_);
  vep_info.hom_lof = vepCount(found_hom_variants, LOF_VEP_FIELD_, LOF_HC_VALUE_);

  vep_info.all_high_effect = vepCount(found_all_unphased, IMPACT_VEP_FIELD_, IMPACT_HIGH_VALUE_);
  vep_info.hom_high_effect = vepCount(found_hom_variants, IMPACT_VEP_FIELD_, IMPACT_HIGH_VALUE_);

  vep_info.all_moderate_effect = vepCount(found_all_unphased, IMPACT_VEP_FIELD_, IMPACT_MODERATE_VALUE_);
  vep_info.hom_moderate_effect = vepCount(found_hom_variants, IMPACT_VEP_FIELD_, IMPACT_MODERATE_VALUE_);

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

        ExecEnv::log().error("GeneVariants::processSummaryStatistics; found unexpected non-zero upper tail_ statistic: {}, population: {}, ",
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

        ExecEnv::log().error("GeneVariants::processSummaryStatistics; found unexpected non-zero lower tail_ statistic: {}, population: {}",
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
