//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to a base count.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::CountFilter::implementFilter(const Variant& variant) const {

  std::shared_ptr<const CountEvidence> count_evidence_ptr = std::dynamic_pointer_cast<const CountEvidence>(variant.evidence());

  if (count_evidence_ptr) {

    return count_evidence_ptr->DPCount() >= minimum_count_;

  } else {

    ExecEnv::log().info("CountFilter; variant does not have base count evidence");

  }

  return true;

}


std::string kgl::CountFilter::filterName() const {

  std::stringstream ss;
  ss << "Variants with minimum base count:" << minimum_count_;
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


