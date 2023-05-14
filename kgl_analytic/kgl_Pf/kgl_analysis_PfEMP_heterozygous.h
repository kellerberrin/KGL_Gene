//
// Created by kellerberrin on 18/04/23.
//

#ifndef KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
#define KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H


#include "kgl_variant_db_population.h"
#include "kgl_pf7_fws_parser.h"
#include "kgl_pf7_sample_parser.h"
#include "kgl_Pf7_physical_distance.h"


#include <memory>



namespace kellerberrin::genome {   //  organization::project level namespace

struct VariantAnalysisType{

  size_t total_variants_{0};
  size_t snp_count_{0};
  size_t indel_count_{0};
  size_t homozygous_count_{0};
  size_t heterozygous_count_{0};
  size_t single_variant_{0};

};

// Indexed by contig.
using VariantAnalysisContigMap = std::map<std::string, VariantAnalysisType>;
class VariantAnalysisObj {

public:

  VariantAnalysisObj(GenomeId_t genome, const Pf7SampleRecord& sample_record, double fws_value)
    : genome_(genome), sample_record_(sample_record), fws_value_(fws_value) {}
  ~VariantAnalysisObj() = default;

  [[nodiscard]] VariantAnalysisContigMap& getMap() { return analysis_map_; }
  [[nodiscard]] const VariantAnalysisContigMap& getConstMap() const { return analysis_map_; }
  [[nodiscard]] double getFWS() const { return fws_value_; }
  [[nodiscard]] double getFIS() const { return fis_value_; }
  [[nodiscard]] const GenomeId_t& getGenome() const { return genome_; }
  [[nodiscard]] const std::string& getCity() const { return sample_record_.location1_; }
  [[nodiscard]] const std::string& getCountry() const { return sample_record_.country_; }
  [[nodiscard]] const std::string& getStudy() const { return sample_record_.study_; }
  [[nodiscard]] const std::string& getYear() const { return sample_record_.year_; }

  void setFIS(double fis_value) { fis_value_ = fis_value; }

private:

  GenomeId_t genome_;
  Pf7SampleRecord sample_record_;
  double fws_value_;
  double fis_value_{0.0};  // Wrights inbreeding coefficient
  VariantAnalysisContigMap analysis_map_;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Location Specific

struct LocationSummary {

  std::string location_;
  LocationType location_type_;
  std::string city_;
  std::string country_;
  std::string region_;
  double  radius_km_;
  size_t radii_samples_;
  size_t radii_samples_OK_;
  std::map<std::string, size_t> studies_;
  double monoclonal_Fst_;
  double hom_het_ratio_;
  size_t total_variants_;
  double variant_rate_;
  size_t single_variant_;    // Heterozygous (A,a)
  size_t homozygous_count_;  // Homozygous (A,A)
  size_t heterozygous_count_; // heterozygous (a, b)
  size_t snp_count_;
  size_t indel_count_;

};

// Key is location.
using LocationSummaryMap = std::map<std::string, LocationSummary>;
// The top level map key is Genome (sample id), the second level is Contig id.
using VariantAnalysisMap = std::map<std::string, VariantAnalysisObj>;
class HeteroHomoZygous {

public:

  HeteroHomoZygous() = default;
  ~HeteroHomoZygous() = default;

  void analyzeVariantPopulation(const std::shared_ptr<const PopulationDB> &gene_population_ptr,
                                const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr,
                                const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr);
  void write_results(const std::string& file_name);
  // If sample vector is empty then summarize all available genomes (samples).
  [[nodiscard]] VariantAnalysisType aggregateResults(const std::vector<GenomeId_t>& sample_vector) const;

  [[nodiscard]] LocationSummaryMap location_summary(const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr,
                                                    const std::shared_ptr<const Pf7SampleLocation>& Pf7_physical_distance_ptr,
                                                    double radius_km,
                                                    const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr) const;

  void write_location_results(const std::string& file_name, const LocationSummaryMap& location_summary) const;

  void UpdateSampleLocation(const LocationSummaryMap& location_summary);

private:

  VariantAnalysisMap variant_analysis_map_;
  constexpr static const char CSV_DELIMITER_ = ',';
  constexpr static const size_t MINIMUM_LOCATION_SAMPLES_ = 20;

};


} // Namespace

#endif //KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
