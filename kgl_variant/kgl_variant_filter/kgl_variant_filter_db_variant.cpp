//
// Created by kellerberrin on 16/10/17.
//

#include "kgl_variant_filter_db_variant.h"


namespace kgl = kellerberrin::genome;


/// Filter variants to a minimum Ref+Alt base count.
bool kgl::RefAltCountFilter::implementFilter(const Variant& variant) const {

  if (auto const& format_data = variant.evidence().formatData()) {

    auto const& fd = *format_data;
    return (fd->refCount() + fd->altCount()) >= minimum_count_;

  }

  ExecEnv::log().info("RefAltCountFilter; variant does not have Ref+Alt base count evidence");

  return true;

}


/// Filter variants to a minimum DP base count.
bool kgl::DPCountFilter::implementFilter(const Variant& variant) const {

  if (auto const& format_data = variant.evidence().formatData()) {

    return (*format_data)->DPCount() >= minimum_count_;

  }

  ExecEnv::log().info("DPCountFilter; variant does not have base count evidence");

  return true;

}


/// Filter to indels that are not mod3 in size (frameshift).
bool kgl::FrameShiftFilter::implementFilter(const Variant& variant) const {

  if (variant.isSNP()) {

    return false;

  }

  return (variant.modifyInterval().second.size() % 3) != 0;

}
