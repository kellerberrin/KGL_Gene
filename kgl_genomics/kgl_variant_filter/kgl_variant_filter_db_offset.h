//
// Created by kellerberrin on 11/08/23.
//

#ifndef KGL_VARIANT_FILTER_DB_OFFSET_H
#define KGL_VARIANT_FILTER_DB_OFFSET_H


#include "kgl_variant_filter_type.h"
#include "kgl_variant_db_genome.h"
#include "kel_utility.h"



namespace kellerberrin::genome {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter homozygous variants. For Example if there are 6 variants defined at the offset consisting of
// two homozygous pairs and 2 singleton (heterozygous) variants. Only the homozygous pairs will be returned.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HomozygousFilter: public FilterOffsets {

public:

  HomozygousFilter() { filterName("Homozygous"); }
  ~HomozygousFilter() override = default;

  // Implemented at offset level
  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<HomozygousFilter>(); }


private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter heterozygous variants. For Example if there are 6 variants defined at the offset consisting of
// two homozygous pairs and 2 singleton heterozygous variants. Only the heterozygous (singleton) variants will be returned.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HeterozygousFilter: public FilterOffsets {

public:

  HeterozygousFilter() { filterName("Heterozygous"); }
  ~HeterozygousFilter() override = default;

  // Implemented at offset level
  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<HeterozygousFilter>(); }


private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Ensure max 2 variants per offset.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DiploidFilter: public FilterOffsets {

public:

  DiploidFilter() { filterName("DiploidFilter"); }
  ~DiploidFilter() override = default;

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<DiploidFilter>(); }

private:

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unique variants disregarding phase. For example, if homozygous then filter to a single variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UniqueUnphasedFilter: public FilterOffsets {

public:

  UniqueUnphasedFilter() { filterName("UniqueUnphasedFilter"); }
  ~UniqueUnphasedFilter() override = default;

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<UniqueUnphasedFilter>(*this); }


private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unique variants including phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UniquePhasedFilter: public FilterOffsets {

public:

  UniquePhasedFilter() { filterName("UniquePhasedFilter"); }
  UniquePhasedFilter(const UniquePhasedFilter&) = default;

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<UniquePhasedFilter>(*this); }

private:


};




} // End namespace.


#endif //KGL_VARIANT_FILTER_DB_OFFSET_H
