//
// Created by kellerberrin on 2/03/18.
//

#ifndef KGL_PLOIDY_ANALYSIS_H
#define KGL_PLOIDY_ANALYSIS_H



#include "kgl_variant_db.h"
#include "kgl_variant_phasing_statistics.h"




struct PloidyData {

  size_t homozygous_{0};
  size_t hq_homozygous_{0};
  size_t heterozygous_{0};
  size_t hq_heterozygous_{0};
  double sum_hq_ratio_{0.0};

};

using PloidyDataMap = std::map<std::string, PloidyData>;

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Ploidy Analysis
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class PloidyAnalysis : public PhasedPopulation {


public:

  PloidyAnalysis(const std::string& analysis) : PhasedPopulation(analysis) {}
  ~PloidyAnalysis() override = default;

  bool addPloidyRecord(const std::string& genome,
                       bool homozygous,
                       bool hq_homozygous,
                       bool heterozygous,
                       bool hq_heterozygous,
                       double ratio);

  bool writePloidyResults(const std::string& file_name, char delimiter ) const;

  constexpr static const char CSV_DELIMITER_ = ',';

private:

  constexpr static const size_t PLOIDY_ARRAY_SIZE_ = 100;

  PloidyDataMap ploidy_data_map_;
  std::array<PloidyData, PLOIDY_ARRAY_SIZE_> ratio_array_;

  mutable std::mutex ploidy_mutex_;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_PLOIDY_ANALYSIS_H
