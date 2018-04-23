//
// Created by kellerberrin on 22/11/17.
//


#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound delete variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




bool kgl::CompoundDeleteFactory::selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const {

  return variant_ptr->isDelete();

}


std::shared_ptr<kgl::Variant> kgl::CompoundDeleteFactory::createCompoundVariant(const CompoundVariantMap& variant_map) const {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  Phred_t quality = calculateQuality(variant_map);

  // create the variant
  std::shared_ptr<Variant> compound_delete(std::make_shared<CompoundDelete>(variant_map.begin()->second->genomeId(),
                                                                            variant_map.begin()->second->contigId(),
                                                                            VariantSequence::UNPHASED,
                                                                            variant_map.begin()->second->offset(),
                                                                            quality,
                                                                            variant_map));

  return compound_delete;

}
