//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_OFFSET_H
#define KGL_VARIANT_DB_OFFSET_H


#include "kgl_variant.h"

#include <vector>
#include <forward_list>
#include <memory>


namespace kellerberrin::genome {   //  organization level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the interface class for the low-level variant container.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using OffsetVariantArray = std::vector<std::shared_ptr<const Variant>>;

class VariantArray {

public:

  VariantArray() = default;
  ~VariantArray() = default;

  VariantArray(const VariantArray &) = delete;
  VariantArray& operator=(const VariantArray &) = delete;


  [[nodiscard]] OffsetVariantArray getVariantArray() const {

    OffsetVariantArray variant_vector = variant_vector_;
    return variant_vector;

  }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) { variant_vector_.emplace_back(variant_ptr); }

  void clearVariant() { variant_vector_.clear(); }


private:

  OffsetVariantArray variant_vector_;

};



} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
