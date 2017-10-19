///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant.h"
#include "kgl_genome_db.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ReadCountFilter : public VariantFilter {

public:

  explicit ReadCountFilter(NucleotideReadCount_t read_count) : read_count_(read_count) {}
  ~ReadCountFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const ReadCountVariant& variant) const override { return variant.readCount() >= read_count_; }
  std::string filterName() const final;

private:

  NucleotideReadCount_t read_count_;

};


class MutantProportionFilter : public VariantFilter {

public:

  explicit MutantProportionFilter(double proportion) : mutant_proportion_(proportion) {}
  ~MutantProportionFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const ReadCountVariant& variant) const override { return variant.proportion() >= mutant_proportion_; }
  std::string filterName() const final;

private:

  double mutant_proportion_;

};

class InCDSFilter : public VariantFilter {

public:

  explicit InCDSFilter(std::shared_ptr<const GenomeDatabase> genome_db_ptr) : genome_db_ptr_(genome_db_ptr) {}
  ~InCDSFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const ReadCountVariant& variant) const override { return implementFilter(variant); }

private:

  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

  bool implementFilter(const Variant& variant) const;

};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
