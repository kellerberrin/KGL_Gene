//
// Created by kellerberrin on 22/11/17.
//


#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound delete variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




bool kgl::CompoundDeleteFactory::selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const {

  if (variant_ptr->isSingle()) {

    const std::shared_ptr<const SNPVariant> SNP_ptr = std::static_pointer_cast<const SNPVariant>(variant_ptr);
    return ExtendDNA5::isDeletion(SNP_ptr->mutant());

  } else {

    return false;

  }

}




std::shared_ptr<const kgl::Variant>
kgl::CompoundDeleteFactory::createCompoundVariant(const CompoundVariantMap& variant_map) const {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  // create the variant
  std::shared_ptr<Variant> compound_delete(std::make_shared<CompoundDelete>(variant_map.begin()->second->variantSource(),
                                                                            variant_map.begin()->second->contig(),
                                                                            variant_map.begin()->second->contigOffset(),
                                                                            variant_map));
  // define its coding sequence.
  compound_delete->defineCoding(variant_map.begin()->second->codingSequences().getFirst());

  return compound_delete;

}
