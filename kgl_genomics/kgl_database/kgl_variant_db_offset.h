//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_OFFSET_H
#define KGL_VARIANT_DB_OFFSET_H


#include "kgl_variant.h"

#include <vector>
#include <list>
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

  OffsetDB() { variant_vector_.reserve(INITIAL_VECTOR_SIZE_); }
  explicit OffsetDB(const OffsetDBArray& update) { variant_vector_ = update; }
  ~OffsetDB() = default;

  OffsetDB(const OffsetDB &) = delete;
  OffsetDB& operator=(const OffsetDB &) = delete;

  [[nodiscard]] const OffsetDBArray& getVariantArray() const { return variant_vector_; }

  void addVariant(std::shared_ptr<const Variant> variant_ptr) { variant_vector_.emplace_back(std::move(variant_ptr)); }
  std::pair<size_t, size_t> inSituFilter(const VariantFilter &filter);

private:

  OffsetDBArray variant_vector_;
  constexpr static const size_t INITIAL_VECTOR_SIZE_ = 2;

  // General purpose filter.
  std::pair<size_t, size_t> inSituGeneral(const VariantFilter &filter);
  // Filters for unique variants up to phase.
  std::pair<size_t, size_t> inSituUnique(const VariantFilter &filter);
  // Ensure max 2 variants per offset.
  // Note that VCF indels are actually offset by -1 because they always contain
  // a reference to the base that preceeds the indel for verification.
  // This code does NOT attempt to adjust for this offset and will only give
  // guaranteed correct results if the DB has been pre-filtered to only contain SNPs.
  std::pair<size_t, size_t> inSituDiploid();


};



} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
