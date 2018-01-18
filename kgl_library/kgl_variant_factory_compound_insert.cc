//
// Created by kellerberrin on 23/11/17.
//


#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound insert variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::CompoundInsertFactory::selectVariant(const std::shared_ptr<const Variant>& variant_ptr) const {

  return variant_ptr->isInsert();

}



std::shared_ptr<kgl::Variant> kgl::CompoundInsertFactory::createCompoundVariant(const CompoundVariantMap& variant_map) const {

  if (variant_map.empty()) {

    ExecEnv::log().error("Variant map is empty");
    return nullptr;

  }

  Phred_t quality = calculateQuality(variant_map);

  // create the variant
  std::shared_ptr<Variant> compound_insert(std::make_shared<CompoundInsert>(variant_map.begin()->second->sourceGenome(),
                                                                            variant_map.begin()->second->contig(),
                                                                            variant_map.begin()->second->contigOffset(),
                                                                            quality,
                                                                            variant_map));

  return compound_insert;

}


