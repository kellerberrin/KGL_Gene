//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_MUTATION_OFFSET_H
#define KGL_VARIANT_MUTATION_OFFSET_H

#include <map>
#include <memory>
#include "kgl_variant_db.h"
#include "kel_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// The indel offset accounting map records previous indels so that inserts and deletes can be properly aligned.
// The returned offset is relative to the offset of the sequence of interest (generally a gene or other region).
using IndelAccountingMap = std::map<ContigOffset_t, SignedOffset_t>;


class AdjustedSequenceOffset {

public:

  AdjustedSequenceOffset() : original_interval_(0, 0), modified_interval_(0, 0)  {}
  AdjustedSequenceOffset(ContigOffset_t sequence_begin_offset, ContigSize_t sequence_size)
  : original_interval_(sequence_begin_offset, sequence_begin_offset + sequence_size),
    modified_interval_(original_interval_) {

    initialSequence(sequence_begin_offset, sequence_size);

  }

  ~AdjustedSequenceOffset() = default;

  // Important - Call this routine before mutation.
  // Important - calculate indel offset adjustment with (this routine) adjustIndelOffsets() BEFORE calling updateIndelAccounting().
  [[nodiscard]] SignedOffset_t adjustIndelOffsets(ContigOffset_t contig_offset) const;

  // Important - Call this routine after mutation.
  // Important - Call this to update the indel offset accounting structure AFTER the actual indel offset has been calculated with adjustIndelOffsets().
  [[nodiscard]] bool updateIndelAccounting(ContigOffset_t contig_offset, SignedOffset_t sequence_size_modify);

  // Total indel offset - can be used to check the sequence size after mutation is complete.
  [[nodiscard]] SignedOffset_t totalIndelOffset() const;

  // Reset the count.
  void initialSequence(ContigOffset_t sequence_begin_offset, ContigSize_t sequence_size) {

    original_interval_.resize(sequence_begin_offset, sequence_begin_offset + sequence_size);
    modified_interval_.resize(sequence_begin_offset, sequence_begin_offset + sequence_size);
    indel_accounting_map_.clear();

  }

  [[nodiscard]] const OpenRightUnsigned& orginalInterval() const { return original_interval_; }
  [[nodiscard]] const OpenRightUnsigned& modifiedInterval() const { return modified_interval_; }

  [[nodiscard]] bool reconcileIntervalOffset() const;

private:

  OpenRightUnsigned original_interval_;
  OpenRightUnsigned modified_interval_;
  IndelAccountingMap indel_accounting_map_;


};



}   // end namespace



#endif // KGL_VARIANT_MUTATION_OFFSET_H
