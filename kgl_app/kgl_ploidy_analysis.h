//
// Created by kellerberrin on 2/03/18.
//

#ifndef KGL_PLOIDY_ANALYSIS_H
#define KGL_PLOIDY_ANALYSIS_H



#include "kgl_variant_db.h"

struct PloidyData {

  size_t homozygous_{0};
  size_t hq_homozygous_{0};
  size_t heterozygous_{0};
  size_t hq_heterozygous_{0};

};

using PloidyDataMap = std::map<std::string, PloidyData>;

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class PloidyAnalysis : public PopulationVariant {


public:

  PloidyAnalysis(const std::string& analysis) : PopulationVariant(analysis) {}
  ~PloidyAnalysis() = default;

  bool addPloidyRecord(const std::string& genome, bool haploid, bool hq_haploid, bool diploid, bool hq_diploid);

  bool writePloidyResults(const std::string& file_name, char delimiter ) const;

  constexpr static const char CSV_DELIMITER_ = ',';

private:

  PloidyDataMap ploidy_data_map_;
  mutable std::mutex ploidy_mutex_;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_KGL_PLOIDY_ANALYSIS_H
