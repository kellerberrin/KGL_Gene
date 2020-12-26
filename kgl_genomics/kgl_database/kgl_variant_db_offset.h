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
// This is the low-level variant container.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

using OffsetDBArray = std::vector<std::shared_ptr<const Variant>>;

class OffsetDB {

public:

  OffsetDB() = default;
  ~OffsetDB() = default;

  OffsetDB(const OffsetDB &) = delete;
  OffsetDB& operator=(const OffsetDB &) = delete;

  [[nodiscard]] const OffsetDBArray& getVariantArray() const { return variant_vector_; }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) { variant_vector_.emplace_back(std::move(variant_ptr)); }

  void clearVariant() { variant_vector_.clear(); }


private:

  OffsetDBArray variant_vector_;

};



} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
