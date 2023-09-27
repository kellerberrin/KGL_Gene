//
// Created by kellerberrin on 21/02/18.
//

#ifndef KGL_RNA_SEARCH_H
#define KGL_RNA_SEARCH_H


#include <memory>
#include <fstream>
#include "kgl_sequence_distance.h"
#include "kgl_sequence_compare.h"
#include "kgl_variant_db_population.h"
#include "kgl_variant_filter_db_variant.h"
#include "kgl_io_gff_fasta.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class RNAAnalysisResults {

public:

  CompareScore_t score_;
  ContigOffset_t rna_offset_;
  StrandSense rna_sense_;
  std::string comparison_;
  StrandSense target_sense_;
  ContigOffset_t target_offset_;
  std::string rna_sequence;

};


using RNASearchResults = std::multimap<CompareScore_t, RNAAnalysisResults>;


class RNAAnalysis {

public:

  RNAAnalysis() = default;
  ~RNAAnalysis() = default;

  [[nodiscard]] bool getRNARegions( const ContigId_t& rna_contig,
                                    ContigOffset_t rna_offset,
                                    ContigSize_t rna_region_size,
                                    StrandSense rna_strand,
                                    const ContigId_t& rna_target_contig,
                                    ContigOffset_t rna_target_offset,
                                    ContigSize_t rna_target_size,
                                    StrandSense rna_target_strand,
                                    const std::shared_ptr<const GenomeReference>& genome_db_ptr);

  [[nodiscard]] bool compareRNARegion( ContigSize_t rna_region_comparison_start,
                                       ContigSize_t rna_region_subsize,
                                       ContigSize_t rna_region_comparison_increment,
                                       const std::shared_ptr<const LocalDNASequenceCompare>& dna_compare_metric);

  void showResults(size_t limit);


private:

  DNA5SequenceLinear rna_sequence_;
  DNA5SequenceLinear rna_target_;
  RNASearchResults search_results_;

};




}   // end namespace genome



#endif //KGL_RNA_SEARCH_H
