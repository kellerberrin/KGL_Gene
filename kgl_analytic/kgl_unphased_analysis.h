//
// Created by kellerberrin on 11/08/18.
//


#ifndef KGL_UNPHASED_ANALYSIS_H
#define KGL_UNPHASED_ANALYSIS_H


#include "kgl_variant_db_unphased_population.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//

class AggregateVariantDistribution {

public:


  explicit AggregateVariantDistribution(ContigSize_t interval) : interval_(interval), analysis_genome_("analysis") {}
  ~AggregateVariantDistribution() = default;

  bool variantDistribution(std::shared_ptr<const UnphasedPopulation> unphased_population);

  bool writeDistribution(std::shared_ptr<const GenomeDatabase> genome_db,
                         const std::string& filename,
                         const char delimiter = ',') const;

private:

  ContigSize_t interval_;
  UnphasedGenome analysis_genome_;

  bool writeData(std::shared_ptr<const GenomeDatabase> genome_db, std::ostream& output, const char delimiter) const;
  bool writeHeader(std::ostream& output, const char delimiter) const;


};


using HeteroResultArray = std::array<size_t, 3>;
using HeteroResults = std::map<GenomeId_t, HeteroResultArray>;
class HeterozygousStatistics {

public:

  HeterozygousStatistics() = default;
  ~HeterozygousStatistics() = default;

  bool heterozygousStatistics(std::shared_ptr<const UnphasedPopulation> unphased_population);
  bool writeHeterozygousStatistics(const std::string& file_name, const char delimiter) const;

private:

  const HeteroResults& getMap() const { return hetero_results_; }

  constexpr static const size_t HOMOZYGOUS_{0};
  constexpr static const size_t HETEROZYGOUS_{1};
  constexpr static const size_t SINGLEHETEROZYGOUS_{2};

  HeteroResults hetero_results_;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_UNPHASED_ANALYSIS_H
