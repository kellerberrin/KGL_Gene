//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include <memory>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;



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


