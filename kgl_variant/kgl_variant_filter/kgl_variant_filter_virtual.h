//
// Created by kellerberrin on 5/2/21.
//

#ifndef KGL_VARIANT_FILTER_VIRTUAL_H
#define KGL_VARIANT_FILTER_VIRTUAL_H

#include <string>
#include <memory>

namespace kellerberrin::genome {   //  organization level namespace

// Runtime filter type enum. Filters are processed differently according to type.
enum class FilterBaseType { POPULATION_FILTER, GENOME_FILTER, CONTIG_FILTER, OFFSET_FILTER, VARIANT_FILTER };


/// Abstract base of every filter. Concrete filters are dispatched by the DB layer on filterType().
///
/// State protocol (load-bearing, enforced by the DB layer in kgl_variant_db):
///   - Variant filters are applied concurrently WITHOUT cloning and must be stateless.
///   - Population/Genome/Contig/Offset filters are cloned by the DB layer before each use and
///     may hold per-invocation mutable state.
class BaseFilter {

public:

  BaseFilter() = default;
  virtual ~BaseFilter() = default;

  [[nodiscard]] std::string filterName() const { return filter_name_; }
  void filterName(std::string filter_name) noexcept { filter_name_ = std::move(filter_name); }

  [[nodiscard]] virtual FilterBaseType filterType() const = 0;
  [[nodiscard]] virtual std::shared_ptr<BaseFilter> clone() const = 0;

private:

  std::string filter_name_;

};



} // namespace

#endif // KGL_VARIANT_FILTER_VIRTUAL_H
