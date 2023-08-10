//
// Created by kellerberrin on 10/08/23.
//

#ifndef KGL_VARIANT_FILTER_DIFFERENCE_H
#define KGL_VARIANT_FILTER_DIFFERENCE_H


#include "kgl_variant_filter_type.h"
#include "kgl_variant_db_genome.h"
#include "kel_utility.h"



namespace kellerberrin::genome {   //  organization::project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference between two genomes.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomeDifferenceFilter : public FilterGenomes {

public:

  explicit GenomeDifferenceFilter(const std::shared_ptr<const GenomeDB>& reference_ptr) : reference_ptr_(reference_ptr) {

    filterName("GenomeDifferenceFilter");

  }
  ~GenomeDifferenceFilter() override = default;

  [[nodiscard]] std::unique_ptr<GenomeDB> applyFilter(const GenomeDB& genome) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<GenomeDifferenceFilter>(*this); }

private:

  std::shared_ptr<const GenomeDB> reference_ptr_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference between two contigs.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigDifferenceFilter : public FilterContigs {

public:

  explicit ContigDifferenceFilter(const std::shared_ptr<const ContigDB>& reference_ptr) : reference_ptr_(reference_ptr) {

    filterName("ContigDifferenceFilter");

  }
  ~ContigDifferenceFilter() override = default;

  [[nodiscard]] std::unique_ptr<ContigDB> applyFilter(const ContigDB& contig) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<ContigDifferenceFilter>(*this); }

private:

  std::shared_ptr<const ContigDB> reference_ptr_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns the difference between two offsets.
// Only variants that do not occur in the reference are returned.
// There is no homozygous or heterozygous difference information returned.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OffsetDifferenceFilter: public FilterOffsets {

public:

  OffsetDifferenceFilter(const std::shared_ptr<const OffsetDB>& reference_ptr) : reference_ptr_(reference_ptr) {

    filterName("OffsetDifferenceFilter");

  }
  OffsetDifferenceFilter(const OffsetDifferenceFilter&) = default;

  [[nodiscard]] std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& offset) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<OffsetDifferenceFilter>(*this); }

private:

  std::shared_ptr<const OffsetDB> reference_ptr_;

};



} // End namespace.


#endif //KGL_VARIANT_FILTER_DIFFERENCE_H
