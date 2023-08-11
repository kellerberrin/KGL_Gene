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

  void addVariant(const std::shared_ptr<const Variant>& variant_ptr) { variant_vector_.push_back(variant_ptr); }

  // Return a filtered copy of the offset.
  // Offset is the bottom level of the population structure, so a true/deep copy of the filtered offset is returned
  // unlike the higher levels of Population/Genome/Contig.
  [[nodiscard]] std::unique_ptr<OffsetDB> viewFilter(const BaseFilter &filter) const;
  // Filter this offset.
  // Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
  std::pair<size_t, size_t> selfFilter(const BaseFilter &filter);

private:

  OffsetDBArray variant_vector_;
  constexpr static const size_t INITIAL_VECTOR_SIZE_ = 2;


};



} // namespace


#endif //KGL_KGL_VARIANT_DB_OFFSET_H
