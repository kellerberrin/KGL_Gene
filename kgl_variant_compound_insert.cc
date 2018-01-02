//
// Created by kellerberrin on 24/11/17.
//


#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;


std::string kgl::CompoundInsert::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;
  ss << "+(" << size() << ")";
  return ss.str() +  location(delimiter, output_index);

}

