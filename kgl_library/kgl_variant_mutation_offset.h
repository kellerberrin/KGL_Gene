//
// Created by kellerberrin on 3/01/18.
//

#ifndef KGL_VARIANT_MUTATION_OFFSET_H
#define KGL_VARIANT_MUTATION_OFFSET_H

#include <map>
#include <memory>
#include "kgl_variant.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// The indel offset accounting map records previous indels so that inserts and deletes can be properly aligned.
using IndelAccountingMap = std::map<ContigOffset_t, SignedOffset_t>;


class VariantMutationOffset {

public:

  explicit VariantMutationOffset() = default;
  ~VariantMutationOffset() = default;

  // Important - Call this routine before mutation.
  // Important - calculate indel offset adjustment with (this routine) adjustIndelOffsets() BEFORE calling updateIndelAccounting().
  SignedOffset_t adjustIndelOffsets(ContigOffset_t contig_offset) const;

  // Important - Call this routine after mutation.
  // Important - Call this to update the indel offset accounting structure AFTER the actual indel offset has been calculated with adjustIndelOffsets().
  bool updateIndelAccounting(std::shared_ptr<const Variant> variant_ptr, SignedOffset_t sequence_size_modify);

  // Total indel offset - can be used to check the sequence size after mutation is complete.
  SignedOffset_t totalIndelOffset() const;

  // Reset the count.
  void clearIndelOffset() { indel_accounting_map_.clear(); }

private:

  IndelAccountingMap indel_accounting_map_;


};





}   // end namespace



#endif // KGL_VARIANT_MUTATION_OFFSET_H
