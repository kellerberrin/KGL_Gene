//
// Created by kellerberrin on 5/2/21.
//

#ifndef KGL_VARIANT_FILTER_VIRTUAL_H
#define KGL_VARIANT_FILTER_VIRTUAL_H


#include <map>
#include <memory>
#include <vector>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract_ VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization level namespace


enum class FilterType   { InfoFilter,
                          VepSubStringFilter,
                          InfoGEQIntegerFilter,
                          InfoGEQFloatFilter,
                          InfoSubStringFilter,
                          InfoBooleanFilter,
                          RefAltCountFilter,
                          DPCountFilter,
                          UniqueUnphasedFilter,  // Uniqueness excludes phase info
                          UniquePhasedFilter,   // Uniqueness includes phase info
                          HomozygousFilter,   // 2 identical variants, phase not tested
                          DiploidFilter, // Ensure max two variants per offset.
                          GenomeFilter, // Filter genomes
                          PhaseFilter,
                          PassFilter,
                          SNPFilter,
                          ContigFilter,
                          VariantFilter,
                          OffsetFilter,
                          RegionFilter,
                          TrueFilter,
                          FalseFilter,
                          NotFilter,
                          AndFilter,
                          OrFilter };

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  [[nodiscard]] virtual FilterType filterType() const = 0;

  [[nodiscard]] std::string filterName() const { return filter_name_; }

  void filterName(std::string filter_name) { filter_name_ = std::move(filter_name); }

  [[nodiscard]] virtual std::shared_ptr<VariantFilter> clone() const = 0;

private:

  std::string filter_name_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomeDB;   // Forward decl.

class FilterGenomes : public VariantFilter {

public:

  FilterGenomes() = default;
  ~FilterGenomes() override = default;

  [[nodiscard]] virtual bool applyFilter(const GenomeDB& genome) const = 0;
  [[nodiscard]] FilterType filterType() const override { return FilterType::GenomeFilter; }

private:

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ContigDB;   // Forward decl.

class FilterContigs : public VariantFilter {

public:

  FilterContigs() = default;
  ~FilterContigs() override = default;

  [[nodiscard]] virtual bool applyFilter(const ContigDB& contig) const = 0;
  [[nodiscard]] FilterType filterType() const override { return FilterType::ContigFilter; }

private:

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OffsetDB;   // Forward decl.

class FilterOffsets : public VariantFilter {

public:

  FilterOffsets() = default;
  ~FilterOffsets() override = default;

  [[nodiscard]] virtual bool applyFilter(const OffsetDB& variant) const = 0;
  [[nodiscard]] FilterType filterType() const override { return FilterType::OffsetFilter; }

private:

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Variant;   // Forward decl.

class FilterVariants : public VariantFilter {

public:

  FilterVariants() = default;
  ~FilterVariants() override = default;

  [[nodiscard]] virtual bool applyFilter(const Variant& variant) const = 0;
  [[nodiscard]] FilterType filterType() const override { return FilterType::VariantFilter; }

private:

};









} // namespace

#endif // KGL_VARIANT_FILTER_VIRTUAL_H
