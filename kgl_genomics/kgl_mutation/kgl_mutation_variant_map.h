//
// Created by kellerberrin on 7/09/23.
//

#ifndef KGL_MUTATION_VARIANT_MAP_H
#define KGL_MUTATION_VARIANT_MAP_H


#include "kgl_mutation_variant_filter.h"
#include "kgl_genome_genome.h"
#include "kgl_variant_db_population.h"
#include "kgl_mutation_analysis.h"


namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Object holds a map suitably sorted and modified variants ready to modify alinear dna sequence over the same [start, end)
// (right open interval with a zero offset).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RegionVariantMap {

public:

  RegionVariantMap() : variant_region_(0, 0) {}
  RegionVariantMap(GenomeId_t genome_id,
                   ContigId_t contig_id,
                   ContigOffset_t start,
                   ContigOffset_t end,
                   OffsetVariantMap variant_map) :
      genome_id_(std::move(genome_id)),
      contig_id_(std::move(contig_id)),
      variant_region_(start, end),
      variant_map_(std::move(variant_map)) {}

  RegionVariantMap(const RegionVariantMap &) = default;

  ~RegionVariantMap() = default;

  [[nodiscard]] const GenomeId_t &genomeId() const { return genome_id_; }
  [[nodiscard]] const ContigId_t &contigId() const { return contig_id_; }
  [[nodiscard]] const OpenRightUnsigned &variantRegion() const { return variant_region_; }
  [[nodiscard]] const OffsetVariantMap &variantMap() const { return variant_map_; }

private:

  GenomeId_t genome_id_;
  ContigId_t contig_id_;
  OpenRightUnsigned variant_region_;
  OffsetVariantMap variant_map_;

};

// Simple struct to return the Region Variant map and some statistics
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RegionReturn {

public:

  RegionReturn(std::shared_ptr<const RegionVariantMap> region_variant_ptr,
               size_t variant_count,
               size_t duplicates,
               size_t upstream_deleted) : region_variant_ptr_(std::move(region_variant_ptr)),
                                          variant_count_(variant_count),
                                          duplicates_(duplicates),
                                          upstream_deleted_(upstream_deleted) {}
    RegionReturn(const RegionReturn&) = default;
    ~RegionReturn() = default;


  [[nodiscard]] const std::shared_ptr<const RegionVariantMap>& regionVariant() const { return region_variant_ptr_; }
  [[nodiscard]] size_t variantCount() const { return variant_count_; }
  [[nodiscard]] size_t duplicates() const { return duplicates_; }
  [[nodiscard]] size_t upstreamDeleted() const { return upstream_deleted_; }

private:

  std::shared_ptr<const RegionVariantMap> region_variant_ptr_;
  size_t variant_count_{0};
  size_t duplicates_{0};
  size_t upstream_deleted_{0};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Detailed info on how a variant updates the sequence interval and sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// The enum below has additional information about delete variants.
// DELETED_REGION - When a large delete variant overlaps the entire region of interest and deletes the region.
// PARTIAL_HIGH_DELETE - When a delete variant begins in the region of interest and deletes 'above' the region.
// PARTIAL_LOW_DELETE - When a delete variant begins below in the region of interest and deletes the beginning of the region.
enum class SequenceUpdateResult { NORMAL, ERROR, DELETED_REGION, PARTIAL_HIGH_DELETE, PARTIAL_LOW_DELETE};

class SequenceVariantUpdate {

public:

  SequenceVariantUpdate(const std::shared_ptr<const Variant>& variant_ptr,
                        const OpenRightUnsigned& prior_interval,
                        const OpenRightUnsigned& post_update_interval,
                        const OpenRightUnsigned& updating_interval,
                        SequenceUpdateResult update_result) :
      variant_ptr_(variant_ptr),
      prior_interval_(prior_interval),
      post_update_interval_(post_update_interval),
      updating_interval_(updating_interval),
      update_result_(update_result) {}
  SequenceVariantUpdate(const SequenceVariantUpdate& copy) = default;
  ~SequenceVariantUpdate() = default;

  [[nodiscard]] const std::shared_ptr<const Variant>& variantPtr() const { return variant_ptr_; }
  [[nodiscard]] const OpenRightUnsigned& priorInterval() const { return prior_interval_; }
  [[nodiscard]] const OpenRightUnsigned& postUpdateInterval() const { return post_update_interval_; }
  [[nodiscard]] const OpenRightUnsigned& updatingInterval() const { return updating_interval_; }
  [[nodiscard]] SequenceUpdateResult updateResult() const { return update_result_; }

  // To determine the actual size of a delete on the interval.
  [[nodiscard]] OpenRightUnsigned priorUpdateIntersect() const { return priorInterval().intersection(updatingInterval()); }

  // String detailing audit information.
  [[nodiscard]] std::string toString() const;

private:

  std::shared_ptr<const Variant> variant_ptr_;
  const OpenRightUnsigned prior_interval_;
  const OpenRightUnsigned post_update_interval_;
  const OpenRightUnsigned updating_interval_;
  SequenceUpdateResult update_result_;

};

using IntervalModifyMap = std::map<ContigOffset_t, SequenceVariantUpdate>;


} // namespace



#endif //KGL_MUTATION_VARIANT_MAP_H
