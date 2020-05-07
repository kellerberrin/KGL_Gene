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

// Stores statistical data per interval

class IntervalData{

public:

  IntervalData() : array_variant_count_(ARRAY_VARIANT_COUNT_, 0) {}
  ~IntervalData() = default;

  void addSNPCount(size_t SNP_count) { SNP_count_ += SNP_count; }
  [[nodiscard]] size_t SNPCount() const { return SNP_count_; }

  void addVariantCount(size_t variant_count) { variant_count_ += variant_count; }
  [[nodiscard]] size_t variantCount() const { return variant_count_; }

  void addArrayVariantCount(size_t size) {

    if (size < ARRAY_VARIANT_COUNT_) {

      ++array_variant_count_[size - 1];

    } else {

      ++array_variant_count_[ARRAY_VARIANT_COUNT_ - 1];

    }

  }
  [[nodiscard]] const std::vector<size_t>& arrayVariantCount() const { return array_variant_count_; }
  void offsetDifference(size_t offset_difference) {

    if (offset_difference > max_offset_difference_) {

      max_offset_difference_ = offset_difference;

    }

    ++active_offset_count_;
    sum_offset_difference_ += offset_difference;

  }
  [[nodiscard]] size_t maxOffsetDifference() const { return max_offset_difference_; }
  [[nodiscard]] double meanOffsetDifference() const {

    if (active_offset_count_ == 0) {

      return 0.0;

    } else {

      return static_cast<double>(sum_offset_difference_) / static_cast<double>(active_offset_count_);

    }

  }

private:

  size_t SNP_count_{0};
  size_t variant_count_{0};
  size_t max_offset_difference_{0};
  size_t active_offset_count_{0};
  size_t sum_offset_difference_{0};

  constexpr static const size_t ARRAY_VARIANT_COUNT_ = 5;
  std::vector<size_t> array_variant_count_;

};


// Analysis class for interval analysis
class IntervalAnalysis : public NullAnalysis {

public:

  IntervalAnalysis() = default;
  ~IntervalAnalysis() = default;

  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "INTERVAL"; }
  [[nodiscard]] std::unique_ptr<NullAnalysis> factory() const override { return std::make_unique<IntervalAnalysis>(); }


  // Setup the analytics to process VCF data.
  // This function must be redefined.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const RuntimeParameterMap& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;

  // Perform the genetic analysis per iteration.
  [[nodiscard]] bool iterateAnalysis( std::shared_ptr<const GenomeCollection> reference_genomes,
                                      std::shared_ptr<const UnphasedPopulation> vcf_iterative_dat) override;

  // All VCF data has been presented, finalize analysis and write results.
  [[nodiscard]] bool finalizeAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes) override;


private:

  constexpr static const char* INTERVAL_SIZE_ = "INTERVALSIZE";
  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";

  size_t interval_size_{0};
  std::string output_file_name_;

  // First count is Variants the second is SNPs.
  using IntervalVector = std::vector<IntervalData>;
  using IntervalMap = std::map<ContigId_t, IntervalVector>;
  IntervalMap interval_map_;

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
