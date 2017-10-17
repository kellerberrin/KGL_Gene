///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ReadCountFilter : public VariantFilter {

public:

  explicit ReadCountFilter(NucleotideReadCount_t read_count) : read_count_(read_count) {}
  ~ReadCountFilter() override = default;

  bool applyFilter(const ReadCountVariant& variant) const final { return variant.readCount() >= read_count_; }
  std::string filterName() const final;

private:

  NucleotideReadCount_t read_count_;

};


class MutantProportionFilter : public VariantFilter {

public:

  explicit MutantProportionFilter(double proportion) : mutant_proportion_(proportion) {}
  ~MutantProportionFilter() override = default;

  bool applyFilter(const ReadCountVariant& variant) const final { return variant.proportion() >= mutant_proportion_; }
  std::string filterName() const final;

private:

  double mutant_proportion_;

};

}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
