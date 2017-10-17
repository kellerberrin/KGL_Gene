//
// Created by kellerberrin on 16/10/17.
//

#include <sstream>
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;

std::string kgl::ReadCountFilter::filterName() const {

  std::ostringstream oss;
  oss << "Minimum Read Count Filter >= " << read_count_;
  return oss.str();

}


std::string kgl::MutantProportionFilter::filterName() const {

  std::ostringstream oss;
  oss << "Mutant Minimum Proportion Filter >= " << mutant_proportion_;
  return oss.str();

}
