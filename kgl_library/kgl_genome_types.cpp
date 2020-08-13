//
// Created by kellerberrin on 18/02/18.
//

#include <sstream>
#include "kgl_genome_types.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant offset output convention.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::offsetOutput(kgl::ContigOffset_t offset, kgl::VariantOutputIndex output_base) {

  std::stringstream ss;

  switch(output_base) {

    case VariantOutputIndex::START_0_BASED:
      ss << offset;
      break;

    case VariantOutputIndex::START_1_BASED:
      ss << (offset + 1);
      break;

  }

  return ss.str();

}
