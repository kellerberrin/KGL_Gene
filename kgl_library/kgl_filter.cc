//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;


bool kgl::InfoFilter::applyFilter(const InfoFilterLambda& filter_lambda, const VCFVariant& variant) const {

  std::optional<std::shared_ptr<const DataMemoryBlock>> info_data_opt = variant.evidence().infoData();

  if (not info_data_opt) {

    return missing_default_;

  }

  auto result = info_data_opt.value()->evidenceHeader()->getSubscribedField(field_name_);

  if (result) {

    InfoDataVariant info_data = result.value().getData(*(variant.evidence().infoData().value()));

    return filter_lambda(info_data);

  } else {

    return missing_default_;

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::RefAltCountFilter::implementFilter(const Variant& variant) const {

  std::optional<std::shared_ptr<const FormatData>> count_evidence_opt = variant.evidence().formatData();

  if (count_evidence_opt) {

    return (count_evidence_opt.value()->refCount() + count_evidence_opt.value()->altCount()) >= minimum_count_;

  } else {

    ExecEnv::log().info("RefAltCountFilter; variant does not have Ref+Alt base count evidence");

  }

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::DPCountFilter::implementFilter(const Variant& variant) const {

  std::optional<std::shared_ptr<const FormatData>> count_evidence = variant.evidence().formatData();

  if (count_evidence) {

    return count_evidence.value()->DPCount() >= minimum_count_;

  } else {

    ExecEnv::log().info("DPCountFilter; variant does not have base count evidence");

  }

  return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigFilter::implementFilter(const Variant& variant) const {

    return variant.contigId() == contig_ident_;

}




