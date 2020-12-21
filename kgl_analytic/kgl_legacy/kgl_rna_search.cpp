//
// Created by kellerberrin on 21/02/18.
//

#include <kgl_sequence/kgl_sequence_offset.h>
#include "kgl_rna_search.h"



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

  // Get the contig.
  std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(rna_contig);
  if (not contig_opt) {

    ExecEnv::log().warn("getRNARegions(), Could not find contig: {} in genome database", rna_contig);
    return false;

  }

  // Check offset and size.
  if ((rna_offset + rna_region_size) > contig_opt.value()->sequence().length() or rna_region_size > contig_opt.value()->sequence().length()) {

    ExecEnv::log().warn("getRNARegions(), contig offset: {} and region size: {} too large for contig: {} length: {}",
                        rna_offset, rna_region_size, contig_opt.value()->contigId(), contig_opt.value()->sequence().length());
    return false;

  }

  // Get the reference DNA sequence
  rna_sequence_ = contig_opt.value()->sequence().subSequence(rna_offset, rna_region_size);
  DNA5SequenceCoding stranded_rna_sequence = SequenceOffset::codingSequence(rna_sequence_, rna_strand);
  rna_sequence_ = DNA5SequenceLinear::downConvertToLinear(stranded_rna_sequence);

  // Get the contig.

  std::optional<std::shared_ptr<const ContigReference>> target_contig_opt = genome_db_ptr->getContigSequence(rna_target_contig);
  if (not target_contig_opt) {

    ExecEnv::log().warn("getRNARegions(), Could not find contig: {} in genome database", rna_target_contig);
    return false;

  }

  // Check offset and size.
  if ((rna_target_offset + rna_target_size) > target_contig_opt.value()->sequence().length() or rna_target_size > target_contig_opt.value()->sequence().length()) {

    ExecEnv::log().warn("getRNARegions(), contig offset: {} and region size: {} too large for contig: {} length: {}",
                        rna_target_offset, rna_target_size, target_contig_opt.value()->contigId(), target_contig_opt.value()->sequence().length());
    return false;

  }

  // Get the RNA target sequence
  rna_target_ = target_contig_opt.value()->sequence().subSequence(rna_target_offset, rna_target_size);
  stranded_rna_sequence = SequenceOffset::codingSequence(rna_target_, rna_target_strand);
  rna_target_ = DNA5SequenceLinear::downConvertToLinear(stranded_rna_sequence);

  return true;

}


bool kgl::RNAAnalysis::compareRNARegion(ContigSize_t rna_region_comparison_start,
                                        ContigSize_t rna_region_subsize,
                                        ContigSize_t rna_region_comparison_increment,
                                        const std::shared_ptr<const LocalDNASequenceCompare>& dna_compare_metric) {


  size_t counter = 0;
  const size_t report_frequency = 10;

  for (size_t idx = rna_region_comparison_start;
       idx < (rna_sequence_.length() - rna_region_subsize);
       idx += rna_region_comparison_increment) {

    DNA5SequenceLinear rna_sub_region = rna_sequence_.subSequence(idx, rna_region_subsize);

    std::string compare_str;
    CompareScore_t score = dna_compare_metric->compare(rna_target_, rna_sub_region, compare_str);

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

