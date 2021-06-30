//
// Created by kellerberrin on 5/2/21.
//

#ifndef KGL_VARIANT_FILTER_VIRTUAL_H
#define KGL_VARIANT_FILTER_VIRTUAL_H


#include <map>
#include <memory>
#include <vector>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
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
                          RegionFilter,
                          TrueFilter,
                          FalseFilter,
                          NotFilter,
                          AndFilter,
                          OrFilter };


class Variant;   // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  [[nodiscard]] virtual bool applyFilter(const Variant& variant) const = 0;

  [[nodiscard]] virtual std::shared_ptr<VariantFilter> clone() const = 0;

  [[nodiscard]] virtual FilterType filterType() const = 0;

  [[nodiscard]] std::string filterName() const { return filter_name_; }

  void filterName(std::string filter_name) { filter_name_ = std::move(filter_name); }


private:

  std::string filter_name_;

};




} // namespace

#endif // KGL_VARIANT_FILTER_VIRTUAL_H
