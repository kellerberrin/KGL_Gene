//
// Created by kellerberrin on 27/04/23.
//

#include "kgl_variant_filter.h"
#include "kgl_variant_filter_Pf7.h"
#include "kgl_analysis_PfEMP.h"
#include "kgl_genome_interval.h"



namespace kgl = kellerberrin::genome;



void kgl::PfEMPAnalysis::countGenes(const std::shared_ptr<const PopulationDB>& population_ptr) {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct OverlappingGenes {

    OverlappingGenes() = default;
    ~OverlappingGenes() = default;
    OverlappingGenes(const OverlappingGenes&) = default;

    IntervalSet coding_interection_;
    std::shared_ptr<const GeneIntervalStructure> gene_struct_ptr_a_;
    std::shared_ptr<const GeneIntervalStructure> gene_struct_ptr_b_;
    std::vector<std::shared_ptr<const Variant>> intersection_variants_;

  };
  using GeneIntersectMap = IntervalMultiMap<std::shared_ptr<OverlappingGenes>>;
  using ContigIntersectMap = std::map<ContigId_t, GeneIntersectMap>;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //
  //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct CountGenes{

    CountGenes(const std::shared_ptr<const GenomeReference>& genome_3D7_ptr) {

      coding_variants_ptr_ = std::make_unique<IntervalCodingVariants>(genome_3D7_ptr);
      createIntersectMap();

    }

    bool countGeneFunc(const std::shared_ptr<const Variant>& variant_ptr) {

      auto gene_vector = coding_variants_ptr_->getGeneCoding(*variant_ptr);

      gene_count_[gene_vector.size()]++;
      countIntersection(variant_ptr);

      return true;

    }

    void countIntersection(const std::shared_ptr<const Variant>& variant_ptr) {

      // Check if the contig is defined.
      auto contig_iter = gene_intersect_map_.find(variant_ptr->contigId());
      if (contig_iter == gene_intersect_map_.end()) {

        return;

      }

      auto const& [contig_id, interval_map] = *contig_iter;
      auto const [offset, extent] = variant_ptr->extentOffset();
      OpenRightInterval variant_interval(offset, offset + extent);

      auto overlap_vector = interval_map.findIntersectsIntervals(variant_interval);
      for (auto& overlap_struct_ptr : overlap_vector) {
        if (overlap_struct_ptr->coding_interection_.intersectsInterval(variant_interval)) {

          ++gene_overlap_;
           overlap_struct_ptr->intersection_variants_.push_back(variant_ptr);

        }

      }

    }

    void createIntersectMap() {

      size_t shared_coding_intervals{0};
      for (auto const& [contig, interval_map] : coding_variants_ptr_->getCodingMap()) {

        auto contig_iter =  gene_intersect_map_.find(contig);
        if (contig_iter == gene_intersect_map_.end()) {

          auto [insert_iter, result] = gene_intersect_map_.try_emplace(contig, GeneIntersectMap());
          if (not result) {

            ExecEnv::log().error("PfEMPAnalysis::CountGenes::createIntersectMap; unable to insert contig: {}", contig);
            return;

          }

          contig_iter = insert_iter;

        }
        // We have a contig map.
        auto& [map_contig, gene_intersect_map] = *contig_iter;

        std::vector<std::tuple<OpenRightInterval, std::shared_ptr<const GeneIntervalStructure>, std::shared_ptr<const GeneIntervalStructure>>> intersect_tuple_vector = interval_map.intersect();
        for (auto const& [interval, gene_a_ptr, gene_b_ptr] : intersect_tuple_vector) {

          auto intersect_set = gene_a_ptr->transcriptUnion().intervalSetIntersection(gene_b_ptr->transcriptUnion());
          if (not intersect_set.empty()) {

            shared_coding_intervals += intersect_set.size();
            OverlappingGenes overlapping_genes;
            overlapping_genes.coding_interection_ = intersect_set;
            overlapping_genes.gene_struct_ptr_a_ = gene_a_ptr;
            overlapping_genes.gene_struct_ptr_b_ = gene_b_ptr;
            gene_intersect_map.insert({interval, std::make_shared<OverlappingGenes>(overlapping_genes)});

          }

        }

      }

    }

    void printIntersects() {

      for (auto const& [contig, intersect_map] :  gene_intersect_map_) {

        for (auto const& [interval, overlapping_ptr] : intersect_map) {

          size_t intersect_size{0};
          for (auto const& intersect_interval : overlapping_ptr->coding_interection_) {

            intersect_size += intersect_interval.size();

          }

          ExecEnv::log().info("PfEMPAnalysis::CountGenes::createIntersectMap; Intersect Gene A: {}, Gene B: {}, Intersect Intervals: {} Intersect Size:{}, Shared Variants: {}",
                              overlapping_ptr->gene_struct_ptr_a_->getGene()->id(),
                              overlapping_ptr->gene_struct_ptr_b_->getGene()->id(),
                              overlapping_ptr->coding_interection_.size(),
                              intersect_size,
                              overlapping_ptr->intersection_variants_.size());
          for (auto const& intersect_interval : overlapping_ptr->coding_interection_) {

            ExecEnv::log().info("PfEMPAnalysis::CountGenes::createIntersectMap; Intersect interval: [{}, {}), interval size: {}",
                                intersect_interval.lower(), intersect_interval.upper(), intersect_interval.size());

          }

        }

      }

    }

    void printResults() {

      for (auto const& [gene_count, variant_count] : gene_count_ ) {

        ExecEnv::log().info("PfEMPAnalysis::countGenes::CountGenes; GeneCount: {}, Variant Count: {}", gene_count, variant_count);

      }
      printIntersects();
      ExecEnv::log().info("countgene::CountGenes::countIntersection; Gene Overlap: {}", gene_overlap_);


    }

    size_t gene_overlap_{0};
    ContigIntersectMap gene_intersect_map_;
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

    filtered_population_ptr = filtered_population_ptr->viewFilter(MultiFilterAllCodingVariants(genome_3D7_ptr_));
    ExecEnv::log().info("Coding Population Final Filtered Size Genome count: {}, Variant Count: {}",
                        filtered_population_ptr->getMap().size(),
                        filtered_population_ptr->variantCount());


  }

  // We need to do a deep copy of the filtered population here since the pass QC and FWS P7 filters only do a shallow copy.
  // And when the resultant population pointers go out of scope they will take the shared population structure with them.
  auto deepcopy_population_ptr = filtered_population_ptr->deepCopy();

  // Filtered population should contain all contigs for all genomes.
  deepcopy_population_ptr->squareContigs();

  ExecEnv::log().info("Population Final Filtered Size Genome count: {}, Variant Count: {}",
                      deepcopy_population_ptr->getMap().size(),
                      deepcopy_population_ptr->variantCount());

  return deepcopy_population_ptr;

}

