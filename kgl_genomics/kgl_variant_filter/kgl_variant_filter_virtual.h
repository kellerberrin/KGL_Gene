//
// Created by kellerberrin on 5/2/21.
//

#ifndef KGL_VARIANT_FILTER_VIRTUAL_H
#define KGL_VARIANT_FILTER_VIRTUAL_H

#include <string>
#include <memory>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract_ VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Runtime Filter Type Enum. Filters are processed differently according to type.
enum class FilterBaseType { POPULATION_FILTER, GENOME_FILTER, CONTIG_FILTER, OFFSET_FILTER, VARIANT_FILTER };


namespace kellerberrin::genome {   //  organization level namespace


class BaseFilter {

public:

  BaseFilter() = default;
  virtual ~BaseFilter() = default;

  [[nodiscard]] std::string filterName() const { return filter_name_; }
  void filterName(std::string filter_name) { filter_name_ = std::move(filter_name); }

  [[nodiscard]] virtual FilterBaseType filterType() const = 0;
  [[nodiscard]] virtual std::shared_ptr<BaseFilter> clone() const = 0;

private:

  std::string filter_name_;

};



} // namespace

#endif // KGL_VARIANT_FILTER_VIRTUAL_H
