//
// Created by kellerberrin on 11/08/18.
//


#ifndef KGL_UNPHASED_ANALYSIS_H
#define KGL_UNPHASED_ANALYSIS_H


#include <kgl_variant_db.h>
#include "kgl_variant_db_unphased_population.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Some simple analysis objects that generate statistics from unphased variants.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using IntervalVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;
using IntervalContigMap = std::map<ContigId_t, IntervalVariantMap>;
class AggregateVariantDistribution {

public:


  explicit AggregateVariantDistribution() = default;
  ~AggregateVariantDistribution() = default;

  bool variantDistribution(std::shared_ptr<const UnphasedPopulation> unphased_population);

  bool variantDistribution(std::shared_ptr<const PhasedPopulation> population_ptr);

  bool writeDistribution(std::shared_ptr<const GenomeDatabase> genome_db,
                         size_t interval_size,
                         const std::string& filename,
                         char delimiter = ',') const;

private:

  IntervalContigMap interval_contig_map_;

  bool writeData(std::shared_ptr<const GenomeDatabase> genome_db, size_t interval_size, std::ostream& output, const char delimiter) const;
  bool writeHeader(std::ostream& output, char delimiter) const;
  bool addVariant(std::shared_ptr<const Variant> variant);

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
