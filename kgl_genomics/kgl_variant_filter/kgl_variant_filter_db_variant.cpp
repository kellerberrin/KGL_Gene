//
// Created by kellerberrin on 16/10/17.
//

#include "kgl_variant_filter_db_variant.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::RefAltCountFilter::implementFilter(const Variant& variant) const {

  if (variant.evidence().formatData()) {

    auto const& format_data = *(variant.evidence().formatData().value());
    return (format_data.refCount() + format_data.altCount()) >= minimum_count_;

  } else {

    ExecEnv::log().info("RefAltCountFilter; variant does not have Ref+Alt base count evidence");

  }

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::DPCountFilter::implementFilter(const Variant& variant) const {


  if (variant.evidence().formatData()) {

    auto const& format_data = *(variant.evidence().formatData().value());
    return format_data.DPCount() >= minimum_count_;

  } else {

    ExecEnv::log().info("DPCountFilter; variant does not have base count evidence");

  }

  return true;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Indels that are not mod3 in size.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::FrameShiftFilter::implementFilter(const Variant& variant) const {

  if (variant.isSNP()) {

    return false;

  }

  const auto [variant_type, modify_interval] = variant.modifyInterval();

  return (modify_interval.size() % 3) != 0;

}



