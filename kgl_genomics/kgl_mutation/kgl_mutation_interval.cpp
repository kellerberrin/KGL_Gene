//
// Created by kellerberrin on 28/08/23.
//

#include "kgl_mutation_interval.h"



namespace kgl = kellerberrin::genome;



std::string kgl::SequenceAuditInfo::toString() const {

  SignedOffset_t indel_size = static_cast<SignedOffset_t>(variant_ptr_->alternateSize())
                              - static_cast<SignedOffset_t>(variant_ptr_->referenceSize());

  std::string audit_string =
      " Prior: " + prior_interval_.toString() +
      " Post: " + post_update_interval_.toString() +
      " Update: " + updating_interval_.toString() +
      " IndelSize: " + std::to_string(indel_size) +
      " Variant: " + variant_ptr_->HGVS();

  return audit_string;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::AdjustedSequenceInterval::updateIndelAccounting(ContigOffset_t allele_offset,
                                                          const std::shared_ptr<const Variant>& variant_ptr,
                                                          const OpenRightInterval& prior_interval,
                                                          const OpenRightInterval& post_update_interval,
                                                          const OpenRightInterval& updating_interval) {

  bool result{true};

  auto const [iter, insert_result] = indel_audit_map_.try_emplace(allele_offset, SequenceAuditInfo(variant_ptr,
                                                                                                   prior_interval,
                                                                                                   post_update_interval,
                                                                                                   updating_interval));
  // Check that the interval has the same size as the accounting map.
  if (not insert_result) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; cannot update (duplicate offset) accounting map, offset: {}, variant:{}",
                        allele_offset, variant_ptr->HGVS());
    result = false;

  }

  if (not reconcileIntervalOffset()) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; accounting map size: {} + original interval size: {} does not equal calculated interval size: {}",
                        intervalSizeModification() , orginalInterval().size(), modifiedInterval().size());

    printAudit();

    result = false;

  }

  return result;

}


kgl::SignedOffset_t kgl::AdjustedSequenceInterval::intervalSizeModification() const {

  SignedOffset_t indel_offset_adjust{0};

  for (auto const& [offset, audit_record] : indel_audit_map_) {

    switch (audit_record.variantPtr()->variantType()) {

      case VariantType::SNP:
        break;

      case VariantType::INDEL_DELETE: {

        OpenRightInterval intersect_update = audit_record.priorInterval().intersection(audit_record.updatingInterval());
        indel_offset_adjust -= static_cast<SignedOffset_t>(intersect_update.size());

        }
        break;

      case VariantType::INDEL_INSERT:
        indel_offset_adjust += static_cast<SignedOffset_t>(audit_record.updatingInterval().size());
        break;

    }

  }

  return indel_offset_adjust;

}

bool kgl::AdjustedSequenceInterval::updateOffsetMap(const std::shared_ptr<const Variant>& variant_ptr) {

  VariantType variant_type = variant_ptr->variantType();

  if (variant_type == VariantType::SNP) {

    return true; // No update.

  }

  ContigOffset_t indel_offset{0};
  OpenRightInterval prior_update_interval(modified_interval_);
  OpenRightInterval updating_interval(0, 0);

  if (variant_type == VariantType::INDEL_DELETE) {

    indel_offset = variant_ptr->offset() + variant_ptr->alternateSize(); // Actual offset of delete.
    ContigOffset_t adjusted_delete_offset = indel_offset + intervalSizeModification(); // Adjusted for previous indels.
    size_t delete_size = variant_ptr->referenceSize() - variant_ptr->alternateSize(); // Size of the delete.
    updating_interval.resize(adjusted_delete_offset, adjusted_delete_offset + delete_size); // Adjusted delete interval
    // Check for that the delete occurs on the modified interval
    if (modified_interval_.disjoint(updating_interval)) {

      // Check that the delete was valid in the first place.
      OpenRightInterval original_delete_interval(indel_offset, indel_offset + delete_size);
      if (not orginalInterval().disjoint(original_delete_interval)) {

        ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; delete variant: {}, delete interval: {} disjoint on modified sequence interval: {}",
                            variant_ptr->HGVS(),
                            updating_interval.toString(),
                            modified_interval_.toString());

        printAudit();

        return false;

      } else { // Original delete was disjoint with the original interval.

        if constexpr (DETAILED_INTERVAL_WARNING_) {

          ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; original delete variant: {}, interval: {} disjoint on  original sequence interval: {}",
                              variant_ptr->HGVS(),
                              original_delete_interval.toString(),
                              orginalInterval().toString());

        }

        return true;

      }

    }

    modified_interval_ = modified_interval_.deleteInterval(updating_interval);


  } else if (variant_type == VariantType::INDEL_INSERT) {

    indel_offset = variant_ptr->offset() + variant_ptr->referenceSize(); // Actual offset of the insert.
    ContigOffset_t adjusted_insert_offset = indel_offset + intervalSizeModification(); // Adjusted for previous indels.
    size_t insert_size = variant_ptr->alternateSize() - variant_ptr->referenceSize(); // Size of the insert.
    updating_interval.resize(adjusted_insert_offset, adjusted_insert_offset + insert_size); // Adjusted insert interval
    if (not modified_interval_.containsOffset(adjusted_insert_offset)) {

      // Check that the insert was valid in the first place.
      if (orginalInterval().containsOffset(indel_offset)) {

        ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; insert variant: {}, insert offset: {} not in sequence interval: {}",
                            variant_ptr->HGVS(),
                            adjusted_insert_offset,
                            modified_interval_.toString());

        printAudit();

        return false;

      } else { // The insert variant was not contained in the original interval

        if constexpr (DETAILED_INTERVAL_WARNING_) {

          ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; insert variant: {}, not in original sequence interval: {}",
                              variant_ptr->HGVS(),
                              original_interval_.toString());

        }

        return true;

      }

    }

    modified_interval_ = modified_interval_.insertInterval(updating_interval);

  } else {

    ExecEnv::log().critical("AdjustedSequenceInterval::updateIndelAccounting; unknown variant type: {} cannot continue",
                            variant_ptr->HGVS());

  }

  return updateIndelAccounting(indel_offset, variant_ptr, prior_update_interval, modified_interval_, updating_interval);

}


bool kgl::AdjustedSequenceInterval::reconcileIntervalOffset() const {

  return static_cast<SignedOffset_t>(modified_interval_.size())
         == (intervalSizeModification() + static_cast<SignedOffset_t>(original_interval_.size()));

}


void kgl::AdjustedSequenceInterval::printAudit() {

  std::scoped_lock lock(audit_mutex_);

  if (audit_count_ < MAX_AUDIT_COUNT_) {

    for (auto const& [offset, audit_record] : indel_audit_map_) {

      ExecEnv::log().info("AdjustedSequenceInterval::printAudit; offset: {}, audit record: {}", offset, audit_record.toString());

    }

  }

  ++audit_count_;

}

