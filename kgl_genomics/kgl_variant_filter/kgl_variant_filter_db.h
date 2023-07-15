//
// Created by kellerberrin on 22/04/23.
//

#ifndef KGL_VARIANT_FILTER_DB_H
#define KGL_VARIANT_FILTER_DB_H

#include "kgl_variant_filter_type.h"
#include "kgl_variant_db_genome.h"
#include "kel_utility.h"



namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters a population to the listed genomes (if they exist).
// Note that this is a shallow copy of the original population.
// Use selfFilter() or deepCopy() to create a permanent population view.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomeListFilter : public FilterPopulations {

public:

  explicit GenomeListFilter(const std::vector<GenomeId_t>& genome_list) : genome_set_(genome_list.begin(), genome_list.end()) {}
  explicit GenomeListFilter(std::set<GenomeId_t> genome_set) : genome_set_(std::move(genome_set)) {}
  ~GenomeListFilter() override = default;

  [[nodiscard]] std::unique_ptr<PopulationDB> applyFilter(const PopulationDB& population) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<GenomeListFilter>(genome_set_); }

private:

  std::set<GenomeId_t> genome_set_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ReturnType both variants if and only there are 2 variants at the location that are identical disregarding phase.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class HomozygousFilter: public FilterOffsets {

public:

  HomozygousFilter() { filterName("Homozygous"); }
  ~HomozygousFilter() override = default;

  // Dummy implementation. Implemented at offset level
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
// Unique variants disregarding phase.
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



#endif //KGL_VARIANT_FILTER_DB_H
