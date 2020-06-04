//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter the info data block.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::InfoFilter::implementFilter(const Variant& variant) const {

  InfoDataEvidence info_evidence_opt = variant.evidence().infoData();

  if (info_evidence_opt) {

    InfoDataVariant info_data = info_field_.getData(*info_evidence_opt.value());

    return filter_lambda_(info_data);

  } else {

    return missing_default_;

  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter an integer field in the info data block.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::InfoGEQIntegerFilter::implementFilter(const Variant& variant) const {

  InfoDataEvidence info_evidence_opt = variant.evidence().infoData();

  if (info_evidence_opt) {

    InfoDataVariant info_data = info_field_.getData(*info_evidence_opt.value());

    return filter_lambda_(info_data);

  } else {

    return missing_default_;

  }

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter an integer field in the info data block.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::InfoGEQFloatFilter::implementFilter(const Variant& variant) const {

  InfoDataEvidence info_evidence_opt = variant.evidence().infoData();

  if (info_evidence_opt) {

    InfoDataVariant info_data = info_field_.getData(*info_evidence_opt.value());

    return filter_lambda_(info_data);

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


std::string kgl::RefAltCountFilter::filterName() const {

  std::stringstream ss;
  ss << "Variants with minimum Ref+Alt base count:" << minimum_count_;
  return ss.str();

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


std::string kgl::DPCountFilter::filterName() const {

  std::stringstream ss;
  ss << "Variants with minimum DP base count:" << minimum_count_;
  return ss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::ContigFilter::implementFilter(const Variant& variant) const {

    return variant.contigId() == contig_ident_;

}

std::string kgl::ContigFilter::filterName() const {

  return "Variants in Contig: " + contig_ident_;

}



std::string kgl::RegionFilter::filterName() const {

  std::stringstream ss;
  ss << "Variant in the half-interval [" << start_ << ", " << end_ << ")";
  return ss.str();

}


