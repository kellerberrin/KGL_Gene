//
// Created by kellerberrin on 27/04/23.
//

#include "kgl_variant_filter.h"
#include "kgl_variant_filter_Pf7.h"
#include "kgl_analysis_PfEMP.h"
#include "kgl_genome_interval.h"



namespace kgl = kellerberrin::genome;


void kgl::PfEMPAnalysis::countGenes(const std::shared_ptr<const PopulationDB>& population_ptr) {

  struct CountGenes{

    CountGenes(const std::shared_ptr<const GenomeReference>& genome_3D7_ptr) {

      coding_variants_ptr_ = std::make_unique<IntervalCodingVariants>(genome_3D7_ptr);

    }

    bool countGeneFunc(const std::shared_ptr<const Variant>& variant_ptr) {

      auto gene_vector = coding_variants_ptr_->getGeneCoding(*variant_ptr);

      gene_count_[gene_vector.size()]++;

      return true;

    }

    void printIntervalStats() {

      size_t shared_coding_intervals{0};
      for (auto const& [contig, interval_map] : coding_variants_ptr_->getCodingMap()) {

        std::vector<std::tuple<OpenRightInterval, std::shared_ptr<const GeneIntervalStructure>, std::shared_ptr<const GeneIntervalStructure>>> intersect_tuple_vector = interval_map.intersect();
        for (auto const& [interval, gene_a_ptr, gene_b_ptr] : intersect_tuple_vector) {

          auto intersect_set = gene_a_ptr->transcriptUnion().intervalSetIntersection(gene_b_ptr->transcriptUnion());
          shared_coding_intervals += intersect_set.size();

        }

        ExecEnv::log().info("PfEMPAnalysis::countGenes::CountGenes; Contig: {}, Intervals (Genes): {}, Overlapping Genes: {}, Shared coding intervals: {}",
                            contig, interval_map.size(), interval_map.intersect().size(), shared_coding_intervals);

      }

    }

    void printResults() {

      printIntervalStats();
      for (auto const& [gene_count, variant_count] : gene_count_ ) {

        ExecEnv::log().info("PfEMPAnalysis::countGenes::CountGenes; GeneCount: {}, Variant Count: {}", gene_count, variant_count);

      }

    }


    std::unique_ptr<IntervalCodingVariants> coding_variants_ptr_;
    std::map<size_t, size_t> gene_count_;

  };

  CountGenes GeneCountObj(genome_3D7_ptr_);
  population_ptr->processAll(GeneCountObj, &CountGenes::countGeneFunc);
  GeneCountObj.printResults();

}


// Quality filter the variants using read depth, VQSLOD and other statistics
std::shared_ptr<kgl::PopulationDB> kgl::PfEMPAnalysis::qualityFilter(const std::shared_ptr<const PopulationDB>& unfiltered_population_ptr) {


  size_t unfiltered_count = unfiltered_population_ptr->variantCount();
  ExecEnv::log().info("Unfiltered Population: {}, Genome count: {}, Variant Count: {}",
                      unfiltered_population_ptr->populationId(),
                      unfiltered_population_ptr->getMap().size(),
                      unfiltered_count);


  auto filtered_qc_population_ptr = unfiltered_population_ptr;
  size_t qc_count{0};
  if constexpr (FILTER_QC_ACTIVE_) {
    // Shallow filter only.
    filtered_qc_population_ptr = Pf7_sample_ptr_->filterPassQCGenomes(unfiltered_population_ptr);

    qc_count = filtered_qc_population_ptr->variantCount();
    double qc_filtered = 100.0 * (static_cast<double>(unfiltered_count - qc_count) / static_cast<double>(unfiltered_count));
    ExecEnv::log().info("Filter P7 QC Pass Genome count: {}, Variant Count: {}, Sample Data Count: {}, filtered: {:.2f}%",
                        filtered_qc_population_ptr->getMap().size(),
                        qc_count,
                        Pf7_sample_ptr_->getMap().size(),
                        qc_filtered);

  } else {

    qc_count = filtered_qc_population_ptr->variantCount();

  }

  auto monoclonal_population_ptr = filtered_qc_population_ptr;
  if constexpr (FILTER_FWS_ACTIVE_) {

    // Shallow filter only.
    monoclonal_population_ptr = Pf7_fws_ptr_->viewFilterFWS(FwsFilterType::GREATER_EQUAL, Pf7FwsResource::MONOCLONAL_FWS_THRESHOLD, filtered_qc_population_ptr);

    size_t monoclonal_count = monoclonal_population_ptr->variantCount();
    double mono_filtered = 100.0 * (static_cast<double>(qc_count - monoclonal_count) / static_cast<double>(qc_count));
    ExecEnv::log().info("Filter MonoClonal FWS: {}, Genome Count: {},Variant Count: {}, filtered: {:.2f}%",
                        Pf7FwsResource::MONOCLONAL_FWS_THRESHOLD,
                        monoclonal_population_ptr->getMap().size(),
                        monoclonal_count,
                        mono_filtered);
  }

  // Call bespoke variant filter.
  P7VariantFilter info_field_filter;
  P7VariantFilter::initializeStats();
  auto filtered_population_ptr = monoclonal_population_ptr->viewFilter(info_field_filter);
  P7VariantFilter::printStats();

  // AF Frequency Filter.
  if constexpr (AF_FILTER_ACTIVE_) {

    P7FrequencyFilter af_frequency_filter(FreqInfoField::AF, VARIANT_FREQUENCY_CUTOFF_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(af_frequency_filter);
    P7FrequencyFilter::printStats(VARIANT_FREQUENCY_CUTOFF_);

  }

  // MLEAF Frequency Filter (redundant).
  if constexpr (MLEAF_FILTER_ACTIVE_) {

    P7FrequencyFilter mleaf_frequency_filter(FreqInfoField::MLEAF, VARIANT_FREQUENCY_CUTOFF_);
    P7FrequencyFilter::initializeStats();
    filtered_population_ptr = filtered_population_ptr->viewFilter(mleaf_frequency_filter);
    P7FrequencyFilter::printStats(VARIANT_FREQUENCY_CUTOFF_);

  }

  // Filter for snp only.
  if constexpr (SNP_FILTER_ACTIVE_) {

    filtered_population_ptr = filtered_population_ptr->viewFilter(SNPFilter());

  }

  if constexpr(CODING_FILTER_ACTIVE_) {

    std::shared_ptr<const PopulationDB> multimap_filtered_population_ptr = filtered_population_ptr->viewFilter(MultiFilterAllCodingVariants(genome_3D7_ptr_));
    ExecEnv::log().info("std::multimap Population Final Filtered Size Genome count: {}, Variant Count: {}",
                        multimap_filtered_population_ptr->getMap().size(),
                        multimap_filtered_population_ptr->variantCount());

    countGenes(multimap_filtered_population_ptr);

  }

  // We need to do a deep copy of the filtered population here since the pass QC and FWS P7 filters only do a shallow copy.
  // And when the resultant population pointers go out of scope they will take the shared population structure with them.
  auto deepcopy_population_ptr = filtered_population_ptr->deepCopy();

  countGenes(deepcopy_population_ptr);

  // Filtered population should contain all contigs for all genomes.
  deepcopy_population_ptr->squareContigs();

  ExecEnv::log().info("Population Final Filtered Size Genome count: {}, Variant Count: {}",
                      deepcopy_population_ptr->getMap().size(),
                      deepcopy_population_ptr->variantCount());

  return deepcopy_population_ptr;

}

