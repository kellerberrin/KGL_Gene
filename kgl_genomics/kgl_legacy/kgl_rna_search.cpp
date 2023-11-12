//
// Created by kellerberrin on 21/02/18.
//

#include "kgl_rna_search.h"

#include <iostream>

namespace kgl = kellerberrin::genome;



bool kgl::RNAAnalysis::getRNARegions(const ContigId_t& rna_contig,
                                     ContigOffset_t rna_offset,
                                     ContigSize_t rna_region_size,
                                     StrandSense rna_strand,
                                     const ContigId_t& rna_target_contig,
                                     ContigOffset_t rna_target_offset,
                                     ContigSize_t rna_target_size,
                                     StrandSense rna_target_strand,
                                     const std::shared_ptr<const GenomeReference>& genome_db_ptr) {

  // Get the contig_ref_ptr.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(rna_contig);
  if (not contig_opt) {

    ExecEnv::log().warn("Could not find contig_ref_ptr: {} in genome database", rna_contig);
    return false;

  }
  const std::shared_ptr<const ContigReference>& contig_ptr = contig_opt.value();

  // Check offset and size.
  OpenRightUnsigned rna_interval(rna_offset, rna_offset+rna_region_size);
  if (contig_ptr->sequence().interval().containsInterval(rna_interval)) {

    ExecEnv::log().warn("RNA interval: {} is not contained in contig: {}, contig_ref_ptr interval: {}",
                        rna_interval.toString(), contig_ptr->contigId(), contig_ptr->sequence().interval().toString());
    return false;

  }

  // Get the reference DNA sequence
  auto rna_sequence_opt = contig_ptr->sequence().subSequence(rna_interval);
  if (not rna_sequence_opt) {

    ExecEnv::log().warn("Cannot extract rna sub-sequence: {} from contig_ref_ptr: {} contig interval: {}",
                        rna_interval.toString(), contig_ptr->contigId(), contig_ptr->sequence().interval().toString());
  }
  const DNA5SequenceLinear& rna_sequence = rna_sequence_opt.value();
  DNA5SequenceCoding stranded_rna_sequence = rna_sequence.codingSequence(rna_strand);
  rna_sequence_ = DNA5SequenceLinear::downConvertToLinear(stranded_rna_sequence);

  // Get the contig_ref_ptr.
  std::optional<std::shared_ptr<const ContigReference>> target_contig_opt = genome_db_ptr->getContigSequence(rna_target_contig);
  if (not target_contig_opt) {

    ExecEnv::log().warn("Could not find contig_ref_ptr: {} in genome database", rna_target_contig);
    return false;

  }
  const std::shared_ptr<const ContigReference>& target_contig_ptr = target_contig_opt.value();

  // Check offset and size.
  OpenRightUnsigned rna_target_interval(rna_target_offset, rna_target_offset + rna_target_size);
  if (target_contig_ptr->sequence().interval().containsInterval(rna_target_interval)) {

    ExecEnv::log().warn("RNA target interval: {} is not contained in target contig: {}, target contig_ref_ptr interval: {}",
                        rna_target_interval.toString(),
                        target_contig_ptr->contigId(),
                        target_contig_ptr->sequence().interval().toString());
    return false;

  }


  // Get the RNA target sequence
  auto rna_target_opt = target_contig_ptr->sequence().subSequence(rna_target_interval);
  if (not rna_target_opt) {

    ExecEnv::log().warn("Cannot extract rna sub-sequence: {} from contig: {} contig_ref_ptr interval: {}",
                        rna_target_interval.toString(),
                        target_contig_ptr->contigId(),
                        target_contig_ptr->sequence().interval().toString());
  }
  const DNA5SequenceLinear& rna_target_sequence = rna_target_opt.value();
  auto stranded_target_rna = rna_target_sequence.codingSequence(rna_target_strand);
  rna_target_ = DNA5SequenceLinear::downConvertToLinear(stranded_target_rna);

  return true;

}


bool kgl::RNAAnalysis::compareRNARegion(ContigSize_t rna_region_comparison_start,
                                        ContigSize_t rna_region_subsize,
                                        ContigSize_t rna_region_comparison_increment,
                                        LinearDistanceMetric dna_compare_metric) {


  size_t counter = 0;
  const size_t report_frequency = 10;

  for (size_t idx = rna_region_comparison_start;
       idx < (rna_sequence_.length() - rna_region_subsize);
       idx += rna_region_comparison_increment) {

    OpenRightUnsigned rna_interval(idx, idx+rna_region_subsize);
    auto rna_sub_region_opt = rna_sequence_.subSequence(rna_interval);
    if (not rna_sub_region_opt) {

      ExecEnv::log().warn("Cannot extract rna sub-region {} from rna sequence: {}",
                          rna_interval.toString(), rna_sequence_.interval().toString());
      return false;

    }
    DNA5SequenceLinear& rna_sub_region = rna_sub_region_opt.value();

    std::string compare_str;
    CompareScore_t score = dna_compare_metric(rna_target_, rna_sub_region);

    CompareScore_t index_score = score * -1;

    RNAAnalysisResults analysis_results;

    analysis_results.rna_offset_ = idx;
    analysis_results.score_ = score;
    analysis_results.target_offset_ = 0;
    analysis_results.comparison_ = compare_str;
    analysis_results.rna_sequence = rna_sub_region.getSequenceAsString();

    search_results_.insert(RNASearchResults::value_type(index_score, analysis_results));

    ++counter;

    if (counter % report_frequency == 0) {

      ExecEnv::log().info("compareRNARegion() just compared: {}, total comparisons: {}",
                          rna_sub_region.getSequenceAsString(), counter);

    }

  }

  return true;

}


void kgl::RNAAnalysis::showResults(size_t limit) {

  size_t counter = 0;
  for (auto result : search_results_) {

    std::cout << "Score:" << result.second.score_ << " RNA Offset:" << result.second.rna_offset_
//              << " " << result.second.rna_sequence
//              << '\n' << result.second.comparison_
              << '\n';

    counter++;

    if (counter >= limit) break;

  }

}

