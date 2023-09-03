//
// Created by kellerberrin on 28/08/23.
//

#include "kgl_mutation_interval.h"



namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::SequenceAuditInfo::toString() const {

  SignedOffset_t indel_size = static_cast<SignedOffset_t>(variant_ptr_->alternateSize())
                              - static_cast<SignedOffset_t>(variant_ptr_->referenceSize());

  std::stringstream ss;
  ss << " Prior: " << prior_interval_.toString()
     << " Post: " << post_update_interval_.toString()
     << " Update: " << updating_interval_.toString()
     << " IndelSize: " << std::to_string(indel_size)
     << " Variant: " << variant_ptr_->HGVS();

  return ss.str();

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

  auto const [variant_type, modify_interval] = variant_ptr->modifyInterval();

  // Adjusted for previous indels.
  OpenRightInterval adjusted_modify_interval(modify_interval);
  adjusted_modify_interval.translate(intervalSizeModification());

  // Save the previous interval
  OpenRightInterval prior_interval(modifiedInterval());

  switch(variant_type) {

    case VariantType::SNP:
      modified_interval_ = updateOffsetSNP(adjusted_modify_interval);
      break;

    case VariantType::INDEL_DELETE:
      modified_interval_ = updateOffsetDelete(adjusted_modify_interval);
      break;

    case VariantType::INDEL_INSERT:
      modified_interval_ = updateOffsetInsert(adjusted_modify_interval);
      break;

  }

  return updateIndelAccounting( modify_interval.lower(),
                                variant_ptr,
                                prior_interval,
                                modified_interval_,
                                adjusted_modify_interval);

}

kel::OpenRightInterval kgl::AdjustedSequenceInterval::updateOffsetSNP(const OpenRightInterval&) {

  return OpenRightInterval{modified_interval_};

}


kel::OpenRightInterval kgl::AdjustedSequenceInterval::updateOffsetInsert(const OpenRightInterval &adj_insert_interval) {

  if (not modified_interval_.containsOffset(adj_insert_interval.lower())) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateOffsetInsert; insert offset: {} not in sequence interval: {}",
                        adj_insert_interval.toString(),
                        modified_interval_.toString());

    printAudit();

    // No update
    return modified_interval_;

  }

  return insertInterval(modified_interval_, adj_insert_interval);

}


kel::OpenRightInterval  kgl::AdjustedSequenceInterval::updateOffsetDelete(const OpenRightInterval &adj_delete_interval) {

  // Check for that the delete occurs on the modified interval
  if (modified_interval_.disjoint(adj_delete_interval)) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; delete interval: {} disjoint on modified sequence interval: {}",
                        adj_delete_interval.toString(),
                        modified_interval_.toString());

    printAudit();

    // No update
    return modified_interval_;

  }

  return deleteInterval(modified_interval_, adj_delete_interval);

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

// To insert an interval the lower() parameter of inserted interval must be within the range [lower, upper).
kel::OpenRightInterval kgl::AdjustedSequenceInterval::insertInterval(const OpenRightInterval& target_interval,
                                                                     const OpenRightInterval &insert_interval) {

  if (target_interval.containsOffset(insert_interval.lower())) {

    return { target_interval.lower(), target_interval.upper() + insert_interval.size()};

  }

  return { target_interval.lower(), target_interval.upper()};

}

// For a valid delete the intersection of the delete interval must be non-empty.
kel::OpenRightInterval kgl::AdjustedSequenceInterval::deleteInterval( const OpenRightInterval& target_interval,
                                                                      const OpenRightInterval &delete_interval) {

  // does not delete this interval.
  if (target_interval.disjoint(delete_interval)) {

    return { target_interval.lower(), target_interval.upper()};

  }

  if (target_interval.containsInterval(delete_interval) and delete_interval.lower() >= target_interval.lower()) {

    return { target_interval.lower(), target_interval.upper()-delete_interval.size()};

  }

  // Not contained and not disjoint therefore obtain the intersection interval
  auto intersect_interval = target_interval.intersection(delete_interval);

  // If the delete interval extends beyond the modified interval.
  if (delete_interval.lower() >= target_interval.lower()) {

    return { target_interval.lower(), target_interval.upper() - intersect_interval.size()};

  }

    // Else the delete interval is below the modified interval.
  size_t upper_adjust = (target_interval.lower() - delete_interval.lower()) + intersect_interval.size();
  return { delete_interval.lower(), target_interval.upper() - upper_adjust};

}

