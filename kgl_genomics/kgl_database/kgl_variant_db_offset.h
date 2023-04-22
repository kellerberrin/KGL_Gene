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
  ~OffsetDB() = default;

  OffsetDB(const OffsetDB &) = delete;
  OffsetDB& operator=(const OffsetDB &) = delete;

  [[nodiscard]] const OffsetDBArray& getVariantArray() const { return variant_vector_; }
  void setVariantArray(const OffsetDBArray& update) { variant_vector_ = update; }

  void addVariant(const std::shared_ptr<const Variant>& variant_ptr) { variant_vector_.push_back(variant_ptr); }

  std::shared_ptr<OffsetDB> offsetFilter(const VariantFilter &filter) const;

  std::pair<size_t, size_t> inSituFilter(const VariantFilter &filter);

  // setIntersection returns an OffsetDB that contains unique variants present in both offsets.
  // The VariantEquality flag determines whether variant phase is used in the equality.
  [[nodiscard]] std::unique_ptr<OffsetDB> setIntersection(const OffsetDB& intersection_offset, VariantEquality variant_equality) const;
  // setComplement returns a OffsetDB that contains variants present in this contig but not present in the complement_contig.
  // This function does NOT guarantee that the returned OffsetDB only has unique variants,
  [[nodiscard]] std::unique_ptr<OffsetDB> setComplement(const OffsetDB& complement_offset, VariantEquality variant_equality) const;
  // setUnion returns a offset that contains the set union of unique variants present in this offset and the union_offset.
  // Specifically, identical variants, defined by the VariantEquality flag as including phase or not, are ignored.
  [[nodiscard]] std::unique_ptr<OffsetDB> setUnion(const OffsetDB& union_offset, VariantEquality variant_equality) const;

private:

  OffsetDBArray variant_vector_;
  constexpr static const size_t INITIAL_VECTOR_SIZE_ = 2;

  // General purpose filter.
  std::pair<size_t, size_t> inSituGeneral(const VariantFilter &filter);
  // Filters for unique variants up to phase.
  std::pair<size_t, size_t> inSituUnique(const VariantFilter &filter);
  // Ensure max 2 variants per offset.
  // Note that VCF indels are generally offset by +1 because they always contain
  // a reference to the base that preceeds the indel for verification.
  std::pair<size_t, size_t> inSituDiploid();
  // If there are 2 identical variants, disregarding phase.
  // Then the variants are retained else they are deleted.
  std::pair<size_t, size_t> inSituHomozygous();


};



} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
