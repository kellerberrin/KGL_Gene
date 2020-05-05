//
// Created by kellerberrin on 5/5/20.
//

#ifndef KGL_ANALYSIS_INTERVAL_H
#define KGL_ANALYSIS_INTERVAL_H

#include "kgl_runtime.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_unphased_population.h"
#include "kgl_analysis_null.h"

namespace kellerberrin::genome {   //  organization::project level namespace



// Stub class for the interval analysis
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

  using IntervalVector = std::vector<std::pair<size_t, size_t>>;
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
