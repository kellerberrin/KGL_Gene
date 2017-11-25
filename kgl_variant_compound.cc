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

  bool result = contigId() == compound_var->contigId()
                and contigOffset() == compound_var->contigOffset()
                and type() == compound_var->type()
                and codingSequenceId() == compound_var->codingSequenceId();

  if (not result) return false;

  if (getMap().size() != compound_var->getMap().size()) return false;

  auto cmp_iterator = compound_var->getMap().begin();
  for (auto iterator = getMap().begin(); iterator != getMap().end(); ++iterator) {

    if(cmp_iterator == compound_var->getMap().end()) return false;

    if (not iterator->second->equivalent(*(cmp_iterator->second))) return false;

    cmp_iterator++;

  }

  return true;

}


