//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A compound variant. A collection of feature aligned and contiguous variants. Insertions and Deletions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::CompoundVariant::equivalent(const Variant& cmp_var) const {

  auto compound_var = dynamic_cast<const CompoundVariant*>(&cmp_var);

  if (compound_var == nullptr) return false;

  for (auto variant : variant_map_) {

    if (not ((*variant.second) == (*this))) return false;

  }

  return contigId() == compound_var->contigId() and contigOffset() == compound_var->contigOffset();

}

