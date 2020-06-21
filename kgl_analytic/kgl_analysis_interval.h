//
// Created by kellerberrin on 5/5/20.
//

#ifndef KGL_ANALYSIS_INTERVAL_H
#define KGL_ANALYSIS_INTERVAL_H

#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"
#include "kgl_analysis_null.h"

#include <array>

namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stores statistical data per interval
class IntervalData{

public:

  IntervalData(ContigId_t contig_id, ContigOffset_t offset, ContigSize_t interval)
  : contig_id_(std::move(contig_id)),
    offset_(offset),
    interval_(interval),
    array_variant_count_(ARRAY_VARIANT_COUNT_, 0) {}
  ~IntervalData() = default;

  [[nodiscard]] ContigId_t contigId() const { return contig_id_; }
  [[nodiscard]] ContigOffset_t offset() const { return offset_; }
  [[nodiscard]] ContigSize_t interval() const { return interval_; }
  // Total SNPs
  void addSNPCount(const size_t SNP_count) { SNP_count_ += SNP_count; }
  [[nodiscard]] size_t SNPCount() const { return SNP_count_; }

  // Total transitions
  void addTransitionCount(const size_t transition_count) { transition_count_ += transition_count; }
  [[nodiscard]] size_t transitionCount() const { return transition_count_; }

  // Variants with vep high and moderate consequences (only Gnomad data)
  void addConsequenceCount(const size_t consequence_count) { consequence_count_ += consequence_count; }
  [[nodiscard]] size_t consequenceCount() const { return consequence_count_; }

  // Total variants
  void addVariantCount(const size_t variant_count) { variant_count_ += variant_count; }
  [[nodiscard]] size_t variantCount() const { return variant_count_; }
  // Number of discrete offsets with variants.
  [[nodiscard]] size_t variantOffsetCount() const { return variant_offset_count_; }

  // No of variants per offset site.
  void addArrayVariantCount(size_t size) ;
  [[nodiscard]] const std::vector<size_t>& arrayVariantCount() const { return array_variant_count_; }

  // Update the offset and variant empty interval.
  void emptyIntervalOffset(const ContigOffset_t& previous_variant_offset, const ContigOffset_t& variant_offset);
  // .first is the offset address (end of the empty zone), .second is the zone size.
  [[nodiscard]] std::pair<ContigOffset_t, ContigSize_t> maxEmptyInterval() const { return max_empty_interval_; }
  // The mean variant empty interval.
  [[nodiscard]] double meanEmptyInterval() const;

private:

  const ContigId_t contig_id_;
  const ContigOffset_t offset_; // offset with the contig.
  const ContigSize_t interval_; // size of this interval.

  size_t SNP_count_{0};
  size_t transition_count_{0};  // The interval Transition count SNPs.
  size_t variant_count_{0};
  size_t consequence_count_{0};
  std::pair<ContigOffset_t , SignedOffset_t> max_empty_interval_{0, 0};
  size_t variant_offset_count_{0};
  size_t sum_empty_interval_{0};

  constexpr static const size_t ARRAY_VARIANT_COUNT_ = 5;
  std::vector<size_t> array_variant_count_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Analysis class for interval analysis
class IntervalAnalysis : public VirtualAnalysis {

public:

  IntervalAnalysis() = default;
  ~IntervalAnalysis() override = default;

  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "INTERVAL"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<IntervalAnalysis>(); }


  // Setup the analytics to process VCF data.
  // This function must be redefined.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const RuntimeParameterMap& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;

  // Perform the genetic analysis per iteration.
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const UnphasedPopulation> vcf_iterative_data) override;

  // Perform the genetic analysis per iteration
  [[nodiscard]] bool iterationAnalysis() override { return true; }

  // All VCF data has been presented, finalize analysis and write results.
  [[nodiscard]] bool finalizeAnalysis() override;


private:

  constexpr static const char* INTERVAL_SIZE_ = "INTERVALSIZE";
  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";

  size_t interval_size_{0};
  std::string output_file_name_;

  // First count is Variants the second is SNPs.
  using IntervalVector = std::vector<IntervalData>;
  using IntervalMap = std::map<ContigId_t, IntervalVector>;
  IntervalMap interval_map_;
  std::shared_ptr<const GenomeReference> genome_;

  constexpr static const char OUTPUT_DELIMITER_ = ',';

  [[nodiscard]] bool getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters);
  void setupIntervalStructure(std::shared_ptr<const GenomeReference> genome);
  [[nodiscard]] bool variantIntervalCount(std::shared_ptr<const UnphasedPopulation> population_ptr);
  [[nodiscard]] bool writeData( std::shared_ptr<const GenomeReference> genome_db, bool display_sequence, std::ostream& output, char delimiter) const;
  [[nodiscard]] bool writeHeader(std::ostream& output, char delimiter, bool display_sequence) const;
  [[nodiscard]] bool writeResults( std::shared_ptr<const GenomeReference> genome_db, bool display_sequence, char delimiter) const;

};



} // namespace


#endif //KGL_ANALYSIS_INTERVAL_H
