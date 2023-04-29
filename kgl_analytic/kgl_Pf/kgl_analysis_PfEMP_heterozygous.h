//
// Created by kellerberrin on 18/04/23.
//

#ifndef KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
#define KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H


#include "kgl_variant_db_population.h"
#include "kgl_pf7_fws_parser.h"

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

using VariantAnalysisContigMap = std::map<std::string, VariantAnalysisType>;
class VariantAnalysisObj {

public:

  VariantAnalysisObj(GenomeId_t genome, double FWS_statistic) : genome_(std::move(genome)), FWS_statistic_(FWS_statistic) {}
  ~VariantAnalysisObj() = default;

  [[nodiscard]] VariantAnalysisContigMap& getMap() { return analysis_map_; }
  [[nodiscard]] double getFWS() const { return FWS_statistic_; }
  [[nodiscard]] const GenomeId_t& getGenome() const { return genome_; }
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
                                const std::shared_ptr<const Pf7FwsResource>& Pf7_fws_ptr);
  void write_results(const std::string& file_name);


private:

  VariantAnalysisMap variant_analysis_map_;
  constexpr static const char CSV_DELIMITER_ = ',';

};


} // Namespace

#endif //KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
