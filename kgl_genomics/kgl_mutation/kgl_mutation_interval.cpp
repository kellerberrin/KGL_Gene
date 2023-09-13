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

// Prints out an update record.
std::string kgl::SequenceVariantUpdate::toString() const {

  SignedOffset_t indel_size = static_cast<SignedOffset_t>(variant_ptr_->alternateSize())
                              - static_cast<SignedOffset_t>(variant_ptr_->referenceSize());


  std::string update_result_text;
  switch(update_result_) {

    case SequenceUpdateResult::NORMAL:
      update_result_text = "Normal";
      break;

    case SequenceUpdateResult::ERROR:
      update_result_text = "Error";
      break;

    case SequenceUpdateResult::DELETED_REGION:
      update_result_text = "Deleted_Region";
      break;

  }

  std::stringstream ss;
  ss << " Update Result: " << update_result_text
     << " Prior: " << prior_interval_.toString()
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


bool kgl::AdjustedSequenceInterval::processVariantMap(const OffsetVariantMap& variant_map) {

  bool result{true};

  for (auto const& [offset, variant_ptr] : variant_map) {

    if (not processVariant(variant_ptr)) {

      result = false;
      ExecEnv::log().warn("AdjustedSequenceInterval::processVariantMap; problem updating with variant: {}, original interval: {} modified interval: {}",
                          variant_ptr->HGVS(), original_interval_.toString(), modified_interval_.toString());

    }

  }

  return result;

}


bool kgl::AdjustedSequenceInterval::updateIndelAccounting(ContigOffset_t allele_offset, SequenceVariantUpdate sequence_update) {

  bool result{true};

  auto const [iter, insert_result] = indel_modify_map_.try_emplace(allele_offset, sequence_update);
  // Check that the interval has the same size as the accounting map.
  if (not insert_result) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; cannot update (duplicate offset) accounting map, offset: {}, variant:{}",
                        allele_offset, sequence_update.variantPtr()->HGVS());
    result = false;

  }

  if (not reconcileIntervalOffset(sequence_update)) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateIndelAccounting; accounting map size: {} + original interval size: {} does not equal calculated interval size: {}",
                        intervalSizeModification() , orginalInterval().size(), modifiedInterval().size());

    printAudit();

    result = false;

  }

  return result;

}


kgl::SignedOffset_t kgl::AdjustedSequenceInterval::intervalSizeModification() const {

  SignedOffset_t indel_size_adjust{0};

  for (auto const& [offset, audit_record] : indel_modify_map_) {

    switch (audit_record.variantPtr()->variantType()) {

      case VariantType::SNP:
        break;

      case VariantType::INDEL_DELETE: {

          OpenRightUnsigned intersect_update = audit_record.priorInterval().intersection(audit_record.updatingInterval());
          indel_size_adjust -= static_cast<SignedOffset_t>(intersect_update.size());

        }
        break;

      case VariantType::INDEL_INSERT:
        indel_size_adjust += static_cast<SignedOffset_t>(audit_record.updatingInterval().size());
        break;

    }

  }

  return indel_size_adjust;

}

kgl::SignedOffset_t kgl::AdjustedSequenceInterval::intervalOffsetModification() const {

  SignedOffset_t indel_offset_adjust{0};

  for (auto const& [offset, audit_record] : indel_modify_map_) {

    switch (audit_record.variantPtr()->variantType()) {

      case VariantType::SNP:
        break;

      case VariantType::INDEL_DELETE:
        indel_offset_adjust -= static_cast<SignedOffset_t>(audit_record.updatingInterval().size());
        break;

      case VariantType::INDEL_INSERT:
        indel_offset_adjust += static_cast<SignedOffset_t>(audit_record.updatingInterval().size());
        break;

    }

  }

  return indel_offset_adjust;

}


bool kgl::AdjustedSequenceInterval::processVariant(const std::shared_ptr<const Variant>& variant_ptr) {

  auto [offset, sequence_audit] = updateInterval(variant_ptr);

  modified_interval_ = sequence_audit.postUpdateInterval();

  return updateIndelAccounting( offset , sequence_audit);

}

std::pair<kgl::ContigOffset_t, kgl::SequenceVariantUpdate>
kgl::AdjustedSequenceInterval::updateInterval(const std::shared_ptr<const Variant>& variant_ptr) const {

  auto const [variant_type, modify_interval] = variant_ptr->modifyInterval();

  // Adjust the modify interval for all previous indel modifications.
  SignedOffset_t previous_offset_modify = intervalOffsetModification();
  auto adjusted_modify_interval = modify_interval.translate(previous_offset_modify);

  // Save the previous interval
  OpenRightUnsigned prior_interval(modifiedInterval());
  OpenRightUnsigned modified_interval(modifiedInterval());

  // Initialize the update result flag.
  SequenceUpdateResult update_result{SequenceUpdateResult::NORMAL};

  // Update according to indel type, no update if an SNP.
  switch(variant_type) {

    case VariantType::SNP: // No interval update for an SNP.
      break;

    case VariantType::INDEL_DELETE: {

        auto updated_interval = updateOffsetDelete(adjusted_modify_interval);
        // No update if an error
        if (updated_interval == modified_interval) {

          adjusted_modify_interval = {adjusted_modify_interval.lower(), adjusted_modify_interval.lower()};
          update_result = SequenceUpdateResult::ERROR;

        } else {

          modified_interval = updated_interval;
          // Check if entire region has been deleted (can happen with small ncRNA genes).
          if (updated_interval.empty()) {

            update_result = SequenceUpdateResult::DELETED_REGION;

          }

        }

      }
      break;

    case VariantType::INDEL_INSERT: {

        auto updated_interval = updateOffsetInsert(adjusted_modify_interval);
        // No update if an error
        if (updated_interval == modified_interval_) {

          // Set the modify interval to zero size.
          adjusted_modify_interval = {adjusted_modify_interval.lower(), adjusted_modify_interval.lower()};
          update_result = SequenceUpdateResult::ERROR;

        } else {

          modified_interval = updated_interval;

        }

      }
      break;

  }

  if (update_result == SequenceUpdateResult::ERROR) {

    ExecEnv::log().warn("AdjustedSequenceInterval::updateInterval; original/modified interval: {}/{}, orig/mod update: {}/{}, variant: {}",
                        original_interval_.toString(),
                        modified_interval_.toString(),
                        modify_interval.toString(),
                        adjusted_modify_interval.toString(),
                        variant_ptr->HGVS());
  }

  SequenceVariantUpdate sequence_variant_update(variant_ptr,
                                                prior_interval,
                                                modified_interval,
                                                adjusted_modify_interval,
                                                update_result);

  return {modify_interval.lower(), sequence_variant_update};

}


kel::OpenRightUnsigned kgl::AdjustedSequenceInterval::updateOffsetInsert(const OpenRightUnsigned &adj_insert_interval) const {

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


kel::OpenRightUnsigned  kgl::AdjustedSequenceInterval::updateOffsetDelete(const OpenRightUnsigned &adj_delete_interval) const {

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


bool kgl::AdjustedSequenceInterval::reconcileIntervalOffset(const SequenceVariantUpdate& sequence_update) const {


  // Explicitly handle the case where the modified interval has been entirely deleted by a large upstream delete variant.
  if (sequence_update.updateResult() == SequenceUpdateResult::DELETED_REGION) {

    return true;

  }

  // Check that the size of the current modified region + all indel size modifications is the same as the original region size.
  SignedOffset_t size_modification = intervalSizeModification();
  return static_cast<SignedOffset_t>(modified_interval_.size())
         == (size_modification + static_cast<SignedOffset_t>(original_interval_.size()));

}

// Used for debugging updates that fail checks.
void kgl::AdjustedSequenceInterval::printAudit() const {

  std::scoped_lock lock(audit_mutex_);

  if (audit_count_ < MAX_AUDIT_COUNT_) {

    for (auto const& [offset, audit_record] : indel_modify_map_) {

      ExecEnv::log().info("AdjustedSequenceInterval::printAudit; offset: {}, audit record: {}", offset, audit_record.toString());

    }

  }

  ++audit_count_;

}

// To insert an interval the lower() parameter of inserted interval must be within the range [lower, upper).
kel::OpenRightUnsigned kgl::AdjustedSequenceInterval::insertInterval(const OpenRightUnsigned& target_interval,
                                                                     const OpenRightUnsigned &insert_interval) {

  if (target_interval.containsOffset(insert_interval.lower())) {

    return { target_interval.lower(), target_interval.upper() + insert_interval.size()};

  }

  return { target_interval.lower(), target_interval.upper()};

}

// The logic of deleting an interval from the target interval is surprisingly complex.
// The logic below should be studied carefully before modification.
kel::OpenRightUnsigned kgl::AdjustedSequenceInterval::deleteInterval(const OpenRightUnsigned& target_interval,
                                                                     const OpenRightUnsigned &delete_interval) {

  // Does not delete this interval.
  if (target_interval.disjoint(delete_interval)) {

    return { target_interval.lower(), target_interval.upper()};

  }

  // Small ncRNA genes can be wholly deleted, check this.
  if (delete_interval.containsInterval(target_interval)) {

    return { delete_interval.lower(), delete_interval.lower() };

  }

  // Check if the delete interval contained in the target interval.
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
  return { delete_interval.lower(), target_interval.upper() - delete_interval.size()};

}

