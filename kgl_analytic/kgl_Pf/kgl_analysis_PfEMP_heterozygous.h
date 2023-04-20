//
// Created by kellerberrin on 18/04/23.
//

#ifndef KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
#define KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H


#include "kgl_variant_db_population.h"

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
// The top level map key is Genome, the second level is Contig id.
using VariantAnalysisMap = std::map<std::string, std::map<std::string, VariantAnalysisType>>;

class HeteroHomoZygous {

public:

  HeteroHomoZygous() = default;
  ~HeteroHomoZygous() = default;

  void analyzeVariantPopulation(const std::shared_ptr<const PopulationDB> &gene_population_ptr);
  void write_results(const std::string& file_name) const;


private:

  VariantAnalysisMap variant_analysis_map_;
  constexpr static const char CSV_DELIMITER_ = ',';

  [[nodiscard]] double getVQSLOD(const std::shared_ptr<const Variant>& variant_ptr);
  constexpr static const char* VQSLOD_FILTER_NAME_{"VQSLOD"};

};


} // Namespace

#endif //KGL_ANALYSIS_PFEMP_HETEROZYGOUS_H
