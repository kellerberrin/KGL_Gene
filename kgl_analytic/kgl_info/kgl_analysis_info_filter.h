//
// Created by kellerberrin on 4/6/20.
//

#ifndef KGL_ANALYSIS_INFO_FILTER_H
#define KGL_ANALYSIS_INFO_FILTER_H


#include "kgl_runtime.h"
#include "kgl_genome_collection.h"
#include "kgl_variant_db_population.h"
#include "kgl_variant_factory_vcf_evidence_analysis.h"

#include "kgl_analysis_null.h"
#include "kgl_filter.h"
#include "kgl_analysis_age.h"

namespace kellerberrin::genome {   //  organization::project level namespace



// Test the Info data filters.
class InfoFilterAnalysis : public VirtualAnalysis {

public:

  InfoFilterAnalysis() = default;
  ~InfoFilterAnalysis() override = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] std::string ident() const override { return "INFO_FILTER"; }
  [[nodiscard]] std::unique_ptr<VirtualAnalysis> factory() const override { return std::make_unique<InfoFilterAnalysis>(); }

  // Setup the analytics to process VCF data.
  // This function must be redefined.
  [[nodiscard]] bool initializeAnalysis( const std::string& work_directory,
                                         const RuntimeParameterMap& named_parameters,
                                         std::shared_ptr<const GenomeCollection> reference_genomes) override;

  // Perform the genetic analysis per VCF file (need not be redefined)
  [[nodiscard]] bool fileReadAnalysis(std::shared_ptr<const DataDB> vcf_population) override;

  // Perform the genetic analysis per iteration (need not be redefined)
  [[nodiscard]] bool iterationAnalysis() override;

  // All VCF data has been presented, finalize analysis and write results (need not be redefined)
  [[nodiscard]] bool finalizeAnalysis() override;

  [[nodiscard]] bool getParameters(const std::string& work_directory, const RuntimeParameterMap& named_parameters);

private:

  std::shared_ptr<PopulationDB> filtered_vcf_population_;
  std::vector<InfoAgeAnalysis> age_analysis_vector_;   // Store for final totals
  AgeSortedMap age_sorted_map_; // final totals.
  std::string output_file_name_;

  constexpr static const char* OUTPUT_FILE_ = "OUTPUTFILE";

  std::shared_ptr<PopulationDB> qualityFilter(std::shared_ptr<const PopulationDB> vcf_population);

  void analyzeField( const std::string& info_field_ident,
                     const std::vector<double>& field_values,
                     std::shared_ptr<const PopulationDB> vcf_population,
                     std::ostream& result_file);

  void analyzeFilteredPopulation( const VariantFilter& filter,
                                  std::shared_ptr<const PopulationDB> vcf_population,
                                  std::ostream& result_file);

  // Analysis by average age.
  void filterByAge(std::shared_ptr<const PopulationDB> vcf_population, std::ostream& result_file);

  bool performAnalysis(std::shared_ptr<const PopulationDB> filtered_population);

  void listAvailableInfoFields(std::shared_ptr<const PopulationDB> vcf_population);


}; // InfoFilterAnalysis








} // namespace


#endif //KGL_ANALYSIS_INFO_FILTER_H
