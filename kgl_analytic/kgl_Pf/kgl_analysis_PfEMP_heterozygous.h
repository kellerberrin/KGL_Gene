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

  VariantAnalysisObj(GenomeId_t genome, double FWS_statistic, std::string city,
                     std::string country, std::string study, std::string year)
  : genome_(std::move(genome)), FWS_statistic_(FWS_statistic), city_(std::move(city)),
    country_(std::move(country)), study_(std::move(study)), year_(std::move(year)) {}
  ~VariantAnalysisObj() = default;

  [[nodiscard]] VariantAnalysisContigMap& getMap() { return analysis_map_; }
  [[nodiscard]] double getFWS() const { return FWS_statistic_; }
  [[nodiscard]] const GenomeId_t& getGenome() const { return genome_; }
  [[nodiscard]] const std::string& getCity() const { return city_; }
  [[nodiscard]] const std::string& getCountry() const { return country_; }
  [[nodiscard]] const std::string& getStudy() const { return study_; }
  [[nodiscard]] const std::string& getYear() const { return year_; }

  bool addContigAnalysis(const ContigId_t& contig, VariantAnalysisType analysis) {

    if (analysis_map_.contains(contig)) {

      ExecEnv::log().error("VariantAnalysisObj::addContigAnalysis; analysis map for genome: {} contains contig: {}", genome_, contig);
      return false;

    }

    analysis_map_[contig] = analysis;
    return true;

  }


private:

  GenomeId_t genome_;
  double FWS_statistic_;
  std::string city_;
  std::string country_;
  std::string study_;
  std::string year_;
  VariantAnalysisContigMap analysis_map_;

};


// The top level map key is Genome, the second level is Contig id.
// using VariantAnalysisMap = std::map<std::string, std::map<std::string, VariantAnalysisType>>;
using VariantAnalysisMap = std::map<std::string, VariantAnalysisObj>;

class HeteroHomoZygous {

public:

  HeteroHomoZygous() = default;
  ~HeteroHomoZygous() = default;

  void analyzeVariantPopulation(const std::shared_ptr<const PopulationDB> &gene_population_ptr,
                                const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr,
                                const std::shared_ptr<const Pf7SampleResource>& Pf7_sample_ptr,
                                const std::shared_ptr<const Pf7SampleLocation>& Pf7_physical_distance_ptr);
  void write_results(const std::string& file_name);


private:

  VariantAnalysisMap variant_analysis_map_;
  constexpr static const char CSV_DELIMITER_ = ',';

};


} // Namespace

#endif //KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
