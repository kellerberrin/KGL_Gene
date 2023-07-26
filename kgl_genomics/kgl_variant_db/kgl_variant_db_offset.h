//
// Created by kellerberrin on 1/7/20.
//

#ifndef KGL_VARIANT_DB_OFFSET_H
#define KGL_VARIANT_DB_OFFSET_H


#include "kgl_variant_db.h"

#include <vector>


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

  OffsetDB(const OffsetDB &) = delete; // Use deepCopy()
  OffsetDB& operator=(const OffsetDB &) = delete; // Use deepCopy()


  [[nodiscard]] const OffsetDBArray& getVariantArray() const { return variant_vector_; }
  void setVariantArray(const OffsetDBArray& update) { variant_vector_ = update; }

  void addVariant(const std::shared_ptr<const Variant>& variant_ptr) { variant_vector_.push_back(variant_ptr); }

  // Return a filtered copy of the offset.
  // Offset is the bottom level of the population structure, so a true/deep copy of the filtered offset is returned
  // unlike the higher levels of Population/Genome/Contig.
  [[nodiscard]] std::unique_ptr<OffsetDB> copyFilter(const BaseFilter &filter) const;
  // Filter this offset.
  // Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> selfFilter(const BaseFilter &filter);
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


};



} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
