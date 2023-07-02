//
// Created by kellerberrin on 22/06/23.
//

#include "kgl_analysis_PfEMP_overlap.h"


namespace kgl = kellerberrin::genome;


kgl::OverlapGenes::OverlapGenes(const std::shared_ptr<const GenomeReference> &genome_ptr) {

  coding_variants_ptr_ = std::make_unique<IntervalCodingVariants>(genome_ptr);
  createIntersectMap();

}


bool kgl::OverlapGenes::countGeneFunc(const std::shared_ptr<const Variant> &variant_ptr) {

  auto gene_vector = coding_variants_ptr_->getGeneCoding(*variant_ptr);
  gene_count_[gene_vector.size()]++;
  countIntersection(variant_ptr);

  return true;

}

void kgl::OverlapGenes::countIntersection(const std::shared_ptr<const Variant> &variant_ptr) {

  // Check if the contig is defined.
  auto contig_iter = gene_intersect_map_.find(variant_ptr->contigId());
  if (contig_iter == gene_intersect_map_.end()) {

    return;

  }

  auto const &[contig_id, interval_map] = *contig_iter;
  auto const [offset, extent] = variant_ptr->extentOffset();
  if (extent == 0) {

    ExecEnv::log().warn("PfEMPAnalysis::OverlapGenes::createIntersectMap; Bad Extent for variant: {}, Cigar: {}, Offset:{}, Extent: {}",
                         variant_ptr->HGVS(), variant_ptr->alternateCigar(), offset, extent);
    return;

  }

  OpenRightInterval variant_interval(offset, offset + extent);

  auto overlap_vector = interval_map.findIntersectsIntervals(variant_interval);
  for (auto &overlap_struct_ptr: overlap_vector) {
    if (overlap_struct_ptr->coding_interection_.intersectsInterval(variant_interval)) {

      if (not overlap_struct_ptr->intersection_variants_.contains(variant_ptr->HGVS())) {

        overlap_struct_ptr->intersection_variants_.insert({variant_ptr->HGVS(), variant_ptr});

      }

    }

  }

}


void kgl::OverlapGenes::createIntersectMap() {

  size_t shared_coding_intervals{0};
  for (auto const &[contig, interval_map]: coding_variants_ptr_->getCodingMap()) {

    auto contig_iter = gene_intersect_map_.find(contig);
    if (contig_iter == gene_intersect_map_.end()) {

      auto [insert_iter, result] = gene_intersect_map_.try_emplace(contig, GeneIntersectMap());
      if (not result) {

        ExecEnv::log().error("PfEMPAnalysis::OverlapGenes::createIntersectMap; unable to insert contig: {}", contig);
        return;

      }

      contig_iter = insert_iter;

    }
    // We have a contig map.
    auto &[map_contig, gene_intersect_map] = *contig_iter;

    std::vector<std::tuple<OpenRightInterval, std::shared_ptr<const GeneIntervalStructure>, std::shared_ptr<const GeneIntervalStructure>>> intersect_tuple_vector = interval_map.intersect();
    for (auto const &[interval, gene_a_ptr, gene_b_ptr]: intersect_tuple_vector) {

      auto intersect_set = gene_a_ptr->transcriptUnion().intervalSetIntersection(gene_b_ptr->transcriptUnion());
      if (not intersect_set.empty()) {

        shared_coding_intervals += intersect_set.size();
        OverlapGenesVariants overlapping_genes;
        overlapping_genes.coding_interection_ = intersect_set;
        overlapping_genes.gene_struct_ptr_a_ = gene_a_ptr;
        overlapping_genes.gene_struct_ptr_b_ = gene_b_ptr;
        gene_intersect_map.insert({interval, std::make_shared<OverlapGenesVariants>(overlapping_genes)});

      }

    }

  }

}

void kgl::OverlapGenes::printIntersects() {

  for (auto const &[contig, intersect_map]: gene_intersect_map_) {

    for (auto const &[interval, overlapping_ptr]: intersect_map) {

      size_t intersect_size{0};
      for (auto const &intersect_interval: overlapping_ptr->coding_interection_) {

        intersect_size += intersect_interval.size();

      }

      ExecEnv::log().info("PfEMPAnalysis::OverlapGenes::createIntersectMap; Intersect Gene A: {}, Gene B: {}, Intersect Intervals: {} Intersect Size:{}, Shared Variants: {}",
                          overlapping_ptr->gene_struct_ptr_a_->getGene()->id(),
                          overlapping_ptr->gene_struct_ptr_b_->getGene()->id(),
                          overlapping_ptr->coding_interection_.size(),
                          intersect_size,
                          overlapping_ptr->intersection_variants_.size());
      for (auto const &intersect_interval: overlapping_ptr->coding_interection_) {

        ExecEnv::log().info("PfEMPAnalysis::OverlapGenes::createIntersectMap; Intersect interval: [{}, {}), interval size: {}",
                            intersect_interval.lower(), intersect_interval.upper(), intersect_interval.size());

      }
      for (auto const &[variant_hash, variant_ptr]: overlapping_ptr->intersection_variants_) {

        auto const [trimmed_reference, trimmed_alternate] = variant_ptr->trimmedSequences();
        auto [offset, extent] = variant_ptr->extentOffset();
        ExecEnv::log().info("PfEMPAnalysis::OverlapGenes::createIntersectMap; Interval Variant: {}, Cigar: {}, Offset: {}, Extent: {}, trimmed Ref: {}, trimmed Alt: {}",
                            variant_hash, variant_ptr->alternateCigar(), offset, extent,
                            trimmed_reference.getSequenceAsString(), trimmed_alternate.getSequenceAsString());

      }

    }

  }

}

void kgl::OverlapGenes::printResults() {

  for (auto const &[gene_count, variant_count]: gene_count_) {

    ExecEnv::log().info("PfEMPAnalysis::countGenes::OverlapGenes; GeneCount: {}, Variant Count: {}", gene_count, variant_count);

  }
  printIntersects();

}


