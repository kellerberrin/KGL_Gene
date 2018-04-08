//
// Created by kellerberrin on 3/04/18.
//

#ifndef KGL_VARIANT_PHASING_STATISTICS_H
#define KGL_VARIANT_PHASING_STATISTICS_H

#include "kgl_variant_factory_vcf_phasing.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object generates phasing statistics. In particular, overlapping variant statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using PhasedSNPVector = std::vector<std::shared_ptr<const Variant>>;
class OffsetPhasingStatistic {

public:

  explicit OffsetPhasingStatistic(ContigOffset_t  offset,
                                  const PhasedSNPVector& phased_snp_vector) : offset_(offset),
                                                                              phased_snp_vector_(phased_snp_vector) {}
  OffsetPhasingStatistic(const OffsetPhasingStatistic&) = default;
  ~OffsetPhasingStatistic() = default;

  ContigOffset_t  offset() const { return offset_; }
  const PhasedSNPVector& snpVector() const { return phased_snp_vector_; }

private:

  ContigOffset_t  offset_;
  PhasedSNPVector phased_snp_vector_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object generates phasing statistics. In particular, overlapping variant statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using ContigStatsMap = std::map<ContigOffset_t, std::shared_ptr<OffsetPhasingStatistic>>;
class ContigPhasingStatistics {

public:

  ContigPhasingStatistics() = default;
  ContigPhasingStatistics(const ContigPhasingStatistics&) = default;
  ~ContigPhasingStatistics() = default;

  bool phasedSNPs(const VCFContigMap& contig);

  const ContigStatsMap& heterozygousSNP() const { return heterozygous_snp_; }
  const ContigStatsMap& homozygousSNP() const { return homozygous_snp_; }
  const ContigStatsMap& singleSNP() const { return single_snp_; }

  bool heterozygousSNP(ContigOffset_t offset, const PhasedSNPVector &phased_snp_vector) { return insertPhasingStatistic(offset, phased_snp_vector, heterozygous_snp_); }
  bool homozygousSNP(ContigOffset_t offset, const PhasedSNPVector &phased_snp_vector) { return insertPhasingStatistic(offset, phased_snp_vector, homozygous_snp_); }
  bool singleSNP(ContigOffset_t  offset, const PhasedSNPVector& phased_snp_vector) { return insertPhasingStatistic(offset, phased_snp_vector, single_snp_); }

private:

  ContigStatsMap heterozygous_snp_;
  ContigStatsMap homozygous_snp_;
  ContigStatsMap single_snp_;

  bool insertPhasingStatistic(ContigOffset_t  offset, const PhasedSNPVector& phased_snp_vector, ContigStatsMap& stats_map);


};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object generates phasing statistics. In particular, overlapping variant statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


using StatContigMap = std::map<ContigId_t, std::shared_ptr<ContigPhasingStatistics>>;
class DiploidPhasingStatistics {

public:

  explicit DiploidPhasingStatistics() = default;
  DiploidPhasingStatistics(const DiploidPhasingStatistics&) = default;
  ~DiploidPhasingStatistics() = default;

  bool phasedSNPs(const VCFGenome& vcf_genome);


  size_t heterozygousSNPCount() const;
  size_t homozygousSNPCount() const;
  size_t singleSNPCount() const;

  const StatContigMap& getMap() const { return contig_map_; }

private:

  StatContigMap contig_map_;


};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This object generates phasing statistics. In particular, overlapping variant statistics.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeStatMap = std::map<GenomeId_t, std::shared_ptr<DiploidPhasingStatistics>>;
class PopulationPhasingStatistics {

public:

  explicit PopulationPhasingStatistics() = default;
  PopulationPhasingStatistics(const PopulationPhasingStatistics&) = default;
  ~PopulationPhasingStatistics() = default;

  bool phasedSNPs(const VCFPopulation& vcf_population);

  size_t phasedSNPCount() const;

  const GenomeStatMap& getMap() const { return genome_map_; }

  bool getPhasing(const GenomeId_t& genome_id,
                  const ContigId_t& contig_id,
                  ContigOffset_t offset,
                  std::shared_ptr<const OffsetPhasingStatistic>& snp_phasing) const;

  void outputPopulation() const;

private:

  GenomeStatMap genome_map_;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_PHASING_STATISTICS_H
