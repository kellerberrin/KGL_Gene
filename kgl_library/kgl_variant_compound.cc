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
                and phaseId() == compound_var->phaseId()
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


// Order variant types.
bool kgl::CompoundVariant::lessThan(const Variant& cmp_var) const {


  if (contigId() < cmp_var.contigId()) {

    return true;

  } else if (contigId() > cmp_var.contigId()) {

    return false;

  } else if (phaseId() < cmp_var.phaseId()) {

    return true;

  } else if (phaseId() > cmp_var.phaseId()) {

    return false;

  } else if (offset() < cmp_var.offset()) {

    return true;

  } else if (offset() > cmp_var.offset()) {

    return false;

  } else if (variantType() < cmp_var.variantType()) {

    return true;

  } else if (variantType() > cmp_var.variantType()) {

    return false;

  }

  auto cmp_compound = dynamic_cast<const CompoundVariant*>(&cmp_var);

  if (not cmp_compound) {

    // Must be a variant type == compound type.
    ExecEnv::log().error("CompoundVariant::lessThan; Expected Compound, got: {}", cmp_var.output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  if (getMap().size() < cmp_compound->getMap().size()) {

    return true;

  } else if (getMap().size() > cmp_compound->getMap().size()) {

    return false;

  }

  auto cmp_iterator = cmp_compound->getMap().begin();
  for (auto variant : getMap()) {

    if(cmp_iterator == cmp_compound->getMap().end()) return false; // should not happen.

    if (variant.second->lessThan(*(cmp_iterator->second))) {

      return true;

    } else if (not variant.second->equivalent(*(cmp_iterator->second))) {

      return false;
    }

    cmp_iterator++;

  }

  return false;

}


void kgl::CompoundVariant::updatePhaseId(PhaseId_t phase_id) {

  protectedPhaseId(phase_id);

  for (auto sub_variant : getMap()) {

    std::shared_ptr<SingleVariant> mutable_sub_variant = std::const_pointer_cast<SingleVariant>(sub_variant.second);
    mutable_sub_variant->updatePhaseId(phase_id);

  }

}


std::string kgl::CompoundVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index) << "\n";

  for (const auto &variant : variant_map_) {

    ss << variant.second->suboutput(delimiter, output_index, detail);

  }

  return ss.str();

}


