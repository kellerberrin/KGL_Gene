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
// ReturnType both variants if and only there are 2 variants at the location that are identical disregarding phase.
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
