//
// Created by kellerberrin on 22/06/23.
//

#ifndef KGL_ANALYSIS_PFEMP_OVERLAP_H
#define KGL_ANALYSIS_PFEMP_OVERLAP_H


#include "kgl_genome_interval.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct OverlapGenesVariants {

  OverlapGenesVariants() = default;
  ~OverlapGenesVariants() = default;
  OverlapGenesVariants(const OverlapGenesVariants &) = default;

  IntervalSet coding_interection_;
  std::shared_ptr<const GeneIntervalStructure> gene_struct_ptr_a_;
  std::shared_ptr<const GeneIntervalStructure> gene_struct_ptr_b_;
  std::map<std::string, std::shared_ptr<const Variant>> intersection_variants_;

};

using GeneIntersectMap = IntervalMultiMap<std::shared_ptr<OverlapGenesVariants>>;
using ContigIntersectMap = std::map<ContigId_t, GeneIntersectMap>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OverlapGenes {

public:

  explicit OverlapGenes(const std::shared_ptr<const GenomeReference> &genome_ptr);
  ~OverlapGenes() = default;

  void printResults();
  bool processPopulation(const std::shared_ptr<const PopulationDB>& population_ptr) { return population_ptr->processAll(*this, &OverlapGenes::countGeneFunc); }

private:

  ContigIntersectMap gene_intersect_map_;
  std::unique_ptr<IntervalCodingVariants> coding_variants_ptr_;
  std::map<size_t, size_t> gene_count_;
  size_t non_canonical_variants_{0};
  size_t canonical_variants_{0};

  bool countGeneFunc(const std::shared_ptr<const Variant> &variant_ptr);
  void countIntersection(const std::shared_ptr<const Variant> &variant_ptr);
  void checkCanonical(const std::shared_ptr<const Variant> &variant_ptr);
  void createIntersectMap();
  void printIntersects();


};


} // namespace


#endif //KGL_ANALYSIS_PFEMP_OVERLAP_H
