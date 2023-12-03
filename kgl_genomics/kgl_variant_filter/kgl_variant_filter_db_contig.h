//
// Created by kellerberrin on 11/08/23.
//

#ifndef KGL_VARIANT_FILTER_DB_CONTIG_H
#define KGL_VARIANT_FILTER_DB_CONTIG_H

#include "kgl_variant_filter_type.h"
#include "kgl_variant_db_genome.h"
#include "kel_interval_unsigned.h"
#include "kel_utility.h"



namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contig Region filter - Uses the half open interval convention [Begin, End).
// Note all contig_ref_ptr offsets begin at 0 (not 1 as is the practice in VCG, GFF3 etc).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigRegionFilter : public FilterContigs {

public:

  explicit ContigRegionFilter(ContigOffset_t start, ContigOffset_t end) : start_(start), end_(end) {

    std::stringstream ss;
    ss << "Contig filter with variants in the half-interval [" << start_ << ", " << end_ << ")";
    filterName(ss.str());

  }
  ~ContigRegionFilter() override = default;

  [[nodiscard]] std::unique_ptr<ContigDB> applyFilter(const ContigDB& contig) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<ContigRegionFilter>(*this); }

private:

  ContigOffset_t start_;
  ContigOffset_t end_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Contig Modify filter - All variants that modify the region [Begin, End).
// Important - this includes upstream indel deletes that extend into the region.
// The offset these deletes is less than start but the delete variants extends into the region.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigModifyFilter : public FilterContigs {

public:

  explicit ContigModifyFilter(ContigOffset_t start, ContigOffset_t end) : start_(start), end_(end) {

    std::stringstream ss;
    ss << "Contig filter with variants that MODIFY DNA in the half-interval [" << start_ << ", " << end_ << ")";
    filterName(ss.str());

  }
  ~ContigModifyFilter() override = default;

  [[nodiscard]] std::unique_ptr<ContigDB> applyFilter(const ContigDB& contig) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<ContigModifyFilter>(*this); }

private:


  ContigOffset_t start_;
  ContigOffset_t end_;

  // A heuristic region [start-margin, start) that looks for indel delete variants upstream of the region [start, end).
  constexpr static const ContigOffset_t UPSTREAM_DOWNSTREAM_MARGIN_{500};

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter out all variants that will be deleted by an upstream Delete variant
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Upstream delete map is indexed by the upper() value of the delete variant regions.
using UpstreamDeleteMap = IntervalUpperMultiMapType<std::shared_ptr<const Variant>>;

class ContigUpstreamFilter : public FilterContigs {

public:

  explicit ContigUpstreamFilter()  { filterName("Contig Upstream Filter"); }
  ~ContigUpstreamFilter() override = default;

  [[nodiscard]] std::unique_ptr<ContigDB> applyFilter(const ContigDB& contig) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<ContigUpstreamFilter>(*this); }

private:

  // Contains all delete variants that are not downstream of another delete variant (are not deleted upstream).
  mutable UpstreamDeleteMap upstream_delete_map_; // Modified in applyFilter().

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns a contig containing all variants in *this contig that match the template contig_ref_ptr.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ContigTemplateFilter : public FilterContigs {

public:

  explicit ContigTemplateFilter(const std::shared_ptr<const ContigDB>& reference_ptr) : reference_ptr_(reference_ptr) {

    filterName("ContigTemplateFilter");

  }
  ~ContigTemplateFilter() override = default;

  [[nodiscard]] std::unique_ptr<ContigDB> applyFilter(const ContigDB& contig) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<ContigTemplateFilter>(*this); }

private:

  std::shared_ptr<const ContigDB> reference_ptr_;

};



} // End namespace.


#endif //KGL_KGL_VARIANT_FILTER_DB_CONTIG_H
