//
// Created by kellerberrin on 28/08/23.
//

#ifndef KGL_MUTATION_INTERVAL_H
#define KGL_MUTATION_INTERVAL_H



#include <map>
#include <memory>
#include "kgl_variant_db.h"
#include "kel_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Detailed info on how a variant updates the sequence interval and sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SequenceAuditInfo {

public:

  SequenceAuditInfo(const std::shared_ptr<const Variant>& variant_ptr,
                    const OpenRightInterval& prior_interval,
                    const OpenRightInterval& post_update_interval,
                    const OpenRightInterval& updating_interval) :
                    variant_ptr_(variant_ptr),
                    prior_interval_(prior_interval),
                    post_update_interval_(post_update_interval),
                    updating_interval_(updating_interval) {}
  SequenceAuditInfo(const SequenceAuditInfo& copy) = default;
  ~SequenceAuditInfo() = default;

  [[nodiscard]] const std::shared_ptr<const Variant>& variantPtr() const { return variant_ptr_; }
  [[nodiscard]] const OpenRightInterval& priorInterval() const { return prior_interval_; }
  [[nodiscard]] const OpenRightInterval& postUpdateInterval() const { return post_update_interval_; }
  [[nodiscard]] const OpenRightInterval& updatingInterval() const { return updating_interval_; }


  // String detailing audit information.
  [[nodiscard]] std::string toString() const;

private:

  std::shared_ptr<const Variant> variant_ptr_;
  const OpenRightInterval prior_interval_;
  const OpenRightInterval post_update_interval_;
  const OpenRightInterval updating_interval_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The indel offset accounting map records previous indels so that inserts and deletes can be properly aligned.
// The returned offset is relative to the offset of the sequence of interest (generally a gene or other region).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using IntervalModifyMap = std::map<ContigOffset_t, SequenceAuditInfo>;

class AdjustedSequenceInterval {

public:

  AdjustedSequenceInterval(ContigOffset_t sequence_begin_offset, ContigSize_t sequence_size)
      : original_interval_(sequence_begin_offset, sequence_begin_offset + sequence_size),
        modified_interval_(original_interval_) {}

  ~AdjustedSequenceInterval() = default;


  [[nodiscard]] const OpenRightInterval& orginalInterval() const { return original_interval_; }
  [[nodiscard]] const OpenRightInterval& modifiedInterval() const { return modified_interval_; }
  [[nodiscard]] const IntervalModifyMap& auditMap() const { return indel_audit_map_; }

  [[nodiscard]] bool updateOffsetMap(const std::shared_ptr<const Variant>& variant_ptr);

private:

  const OpenRightInterval original_interval_;
  OpenRightInterval modified_interval_;
  IntervalModifyMap indel_audit_map_;

  [[nodiscard]] bool reconcileIntervalOffset() const;

  [[nodiscard]] bool updateIndelAccounting(ContigOffset_t contig_offset,
                                           const std::shared_ptr<const Variant>& variant_ptr,
                                           const OpenRightInterval& prior_interval,
                                           const OpenRightInterval& post_update_interval,
                                           const OpenRightInterval& updating_interval);

  // Calculates the sequence size after all mutations.
  [[nodiscard]] SignedOffset_t intervalSizeModification() const;

  [[nodiscard]] OpenRightInterval updateOffsetSNP(const OpenRightInterval &adj_snp_interval);
  [[nodiscard]] OpenRightInterval updateOffsetInsert(const OpenRightInterval &adj_insert_interval);
  [[nodiscard]] OpenRightInterval updateOffsetDelete(const OpenRightInterval &adj_delete_interval);

  // Insert and Delete are used to modify an interval as if modified by the inserted and deleted intervals of indel variants.
  // To insert an interval the lower() parameter of inserted interval must be within the range [lower, upper).
  [[nodiscard]] static OpenRightInterval insertInterval(const OpenRightInterval& target_interval, const OpenRightInterval &insert_interval);
  // For a valid delete the intersection of the delete interval must be non-empty.
  // Note that the delete_interval argument may be modified if it is not fully contained in this interval.
  [[nodiscard]] static OpenRightInterval deleteInterval(const OpenRightInterval& target_interval, const OpenRightInterval &delete_interval);


  // Detailed output for unexpected conditions.
  void printAudit();

  // Interval calculations will probably be multi-threaded, thus we control thread access to the audit output.
  mutable std::mutex audit_mutex_;
  inline static size_t audit_count_{0};
  constexpr static const size_t MAX_AUDIT_COUNT_{10} ; // Limit the number of audits printed to control log size.

};



  }   // end namespace



#endif //KGL_MUTATION_INTERVAL_H
