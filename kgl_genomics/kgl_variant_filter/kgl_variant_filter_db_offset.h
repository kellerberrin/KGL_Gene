//
// Created by kellerberrin on 11/08/23.
//

#ifndef KGL_VARIANT_FILTER_DB_OFFSET_H
#define KGL_VARIANT_FILTER_DB_OFFSET_H


#include "kgl_variant_filter_type.h"
#include "kgl_variant_db_genome.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/// Filter to offsets consisting of exactly two identical (homozygous) variants.
/// All other offset sizes yield an empty offset. For example, an offset with 6 variants consisting
/// of two homozygous pairs and 2 singleton (heterozygous) variants returns an EMPTY offset.
class HomozygousFilter : public FilterOffsets {

public:

  HomozygousFilter() { filterName("Homozygous"); }
  ~HomozygousFilter() override = default;

  // Implemented at offset level
  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<HomozygousFilter>(); }

};


/// Filter to singleton (heterozygous) variants only.
/// For example, if there are 6 variants at the offset consisting of two homozygous pairs and
/// 2 singleton heterozygous variants, only the heterozygous (singleton) variants are returned.
class HeterozygousFilter : public FilterOffsets {

public:

  HeterozygousFilter() { filterName("Heterozygous"); }
  ~HeterozygousFilter() override = default;

  // Implemented at offset level
  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<HeterozygousFilter>(); }

};


/// Ensure max 2 variants per offset.
class DiploidFilter : public FilterOffsets {

public:

  DiploidFilter() { filterName("DiploidFilter"); }
  ~DiploidFilter() override = default;

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<DiploidFilter>(); }

};


/// Unique variants disregarding phase. For example, if homozygous then filter to a single variant.
class UniqueUnphasedFilter : public FilterOffsets {

public:

  UniqueUnphasedFilter() { filterName("UniqueUnphasedFilter"); }
  ~UniqueUnphasedFilter() override = default;

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<UniqueUnphasedFilter>(*this); }

};


/// Unique variants including phase.
class UniquePhasedFilter : public FilterOffsets {

public:

  UniquePhasedFilter() { filterName("UniquePhasedFilter"); }

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<UniquePhasedFilter>(*this); }

};




} // End namespace.


#endif //KGL_VARIANT_FILTER_DB_OFFSET_H
