//
// Created by kellerberrin on 28/08/23.
//

#ifndef KGL_MUTATION_INTERVAL_H
#define KGL_MUTATION_INTERVAL_H



#include <map>
#include <memory>
#include "kgl_variant_db.h"
#include "kgl_mutation_variant_map.h"
#include "kel_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The indel offset accounting map records previous indels so that inserts and deletes can be properly aligned.
// The returned offset is relative to the offset of the sequence of interest (generally a gene or other region).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AdjustedSequenceInterval {

public:

  AdjustedSequenceInterval(const OpenRightUnsigned& original_interval) : original_interval_(original_interval),
                                                                         modified_interval_(original_interval_) {}

  ~AdjustedSequenceInterval() = default;


  [[nodiscard]] const OpenRightUnsigned& orginalInterval() const { return original_interval_; }
  [[nodiscard]] const OpenRightUnsigned& modifiedInterval() const { return modified_interval_; }
  [[nodiscard]] const IntervalModifyMap& indelModifyMap() const { return indel_modify_map_; }

  [[nodiscard]] bool processVariantMap(const OffsetVariantMap& variant_map);

private:

  const OpenRightUnsigned original_interval_;
  OpenRightUnsigned modified_interval_;
  IntervalModifyMap indel_modify_map_;

  [[nodiscard]] bool processVariant(const std::shared_ptr<const Variant>& variant_ptr);

  [[nodiscard]] bool reconcileIntervalOffset(const SequenceVariantUpdate& sequence_update) const;

  [[nodiscard]] std::pair<ContigOffset_t, SequenceVariantUpdate> updateInterval(const std::shared_ptr<const Variant>& variant_ptr) const;

  [[nodiscard]] bool updateIndelAccounting(ContigOffset_t allele_offset, SequenceVariantUpdate sequence_update);

  // Calculates the sequence size after all mutations.
  [[nodiscard]] SignedOffset_t intervalSizeModification() const;
  // Calculates the sequence offset after all mutations. Note this is not the same as size modification.
  [[nodiscard]] SignedOffset_t intervalOffsetModification() const;

  // Returns the updated interval or the unmodified interval if a problem encountered.
  [[nodiscard]] OpenRightUnsigned updateOffsetInsert(const OpenRightUnsigned &adj_insert_interval) const;
  [[nodiscard]] OpenRightUnsigned updateOffsetDelete(const OpenRightUnsigned &adj_delete_interval) const;

  // Insert and Delete are used to modify an interval as if modified by the inserted and deleted intervals of indel variants.
  // To insert an interval the lower() parameter of inserted interval must be within the range [lower, upper).
  [[nodiscard]] static OpenRightUnsigned insertInterval(const OpenRightUnsigned& target_interval, const OpenRightUnsigned &insert_interval);
  // For a valid delete the intersection of the delete interval must be non-empty.
  // Note that the delete_interval argument may be modified if it is not fully contained in this interval.
  [[nodiscard]] static OpenRightUnsigned deleteInterval(const OpenRightUnsigned& target_interval, const OpenRightUnsigned &delete_interval);


  // Detailed output for unexpected conditions.
  void printAudit() const;

  // Interval calculations will probably be multi-threaded, thus we control thread access to the audit output.
  mutable std::mutex audit_mutex_;
  inline static size_t audit_count_{0};
  constexpr static const size_t MAX_AUDIT_COUNT_{10} ; // Limit the number of audits printed to control log size.

};



  }   // end namespace



#endif //KGL_MUTATION_INTERVAL_H
