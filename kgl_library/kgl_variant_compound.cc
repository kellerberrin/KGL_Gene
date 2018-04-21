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
                and offset() == compound_var->offset()
                and variantType() == compound_var->variantType();

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


std::string kgl::CompoundVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << quality() << delimiter;
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index) << "\n";

  for (const auto &variant : variant_map_) {

    ss << variant.second->suboutput(delimiter, output_index, detail);

  }

  return ss.str();

}


