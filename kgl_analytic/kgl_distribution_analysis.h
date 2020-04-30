//
// Created by kellerberrin on 11/08/18.
//


#ifndef KGL_UNPHASED_ANALYSIS_H
#define KGL_UNPHASED_ANALYSIS_H


#include <kgl_variant_db.h>
#include "kgl_variant_db_unphased_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace


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

  [[nodiscard]] bool variantDistribution(std::shared_ptr<const UnphasedPopulation> unphased_population);

  [[nodiscard]] bool variantDistribution(std::shared_ptr<const PhasedPopulation> population_ptr);

  [[nodiscard]] bool writeDistribution( std::shared_ptr<const RuntimeGenomeDatabase> genome_db,
                                        size_t interval_size,
                                        const ContigId_t analysis_contig,
                                        ContigOffset_t start_offset,
                                        ContigOffset_t end_offset,
                                        bool display_sequence,
                                        const std::string& filename,
                                        char delimiter = ',') const;

private:

  IntervalContigMap interval_contig_map_;

  [[nodiscard]] bool writeData( std::shared_ptr<const RuntimeGenomeDatabase> genome_db,
                                size_t interval_size,
                                const ContigId_t analysis_contig,
                                ContigOffset_t start_offset,
                                ContigOffset_t end_offset,
                                bool display_sequence,
                                std::ostream& output,
                                const char delimiter) const;
  [[nodiscard]] bool writeHeader(std::ostream& output, char delimiter, bool display_sequence) const;
  [[nodiscard]] bool addVariant(std::shared_ptr<const Variant> variant);

};


using HeteroResultArray = std::array<size_t, 3>;
using HeteroResults = std::map<GenomeId_t, HeteroResultArray>;
class HeterozygousStatistics {

public:

  HeterozygousStatistics() = default;
  ~HeterozygousStatistics() = default;

  [[nodiscard]] bool heterozygousStatistics(std::shared_ptr<const UnphasedPopulation> unphased_population);
  [[nodiscard]] bool writeHeterozygousStatistics(const std::string& file_name, const char delimiter) const;

private:

  [[nodiscard]]  const HeteroResults& getMap() const { return hetero_results_; }

  constexpr static const size_t HOMOZYGOUS_{0};
  constexpr static const size_t HETEROZYGOUS_{1};
  constexpr static const size_t SINGLEHETEROZYGOUS_{2};

  HeteroResults hetero_results_;

};



}   // end namespace



#endif //KGL_UNPHASED_ANALYSIS_H
