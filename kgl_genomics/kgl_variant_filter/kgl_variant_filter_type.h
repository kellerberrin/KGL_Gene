//
// Created by kellerberrin on 24/04/23.
//

#ifndef KGL_VARIANT_FILTER_TYPE_H
#define KGL_VARIANT_FILTER_TYPE_H

#include "kgl_variant_filter_virtual.h"
#include "kgl_variant_db_population.h"

#include <map>
#include <vector>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Base class for a filter at the PopulationDB level. Filters genomes within a population.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class FilterPopulations : public BaseFilter {

public:

  FilterPopulations() = default;
  ~FilterPopulations() override = default;

  [[nodiscard]] virtual std::unique_ptr<PopulationDB> applyFilter(const PopulationDB& genome) const = 0;
  [[nodiscard]] FilterBaseType filterType() const override { return FilterBaseType::POPULATION_FILTER; }

private:

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Base class for a filter at the GenomeDB level. Filters contigs within a genome.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterGenomes : public BaseFilter {

public:

  FilterGenomes() = default;
  ~FilterGenomes() override = default;

  [[nodiscard]] virtual std::unique_ptr<GenomeDB> applyFilter(const GenomeDB& genome) const = 0;
  [[nodiscard]] FilterBaseType filterType() const override { return FilterBaseType::GENOME_FILTER; }

private:

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Base class for a filter at the ContigDB level. Filters offsets within a contig_ref_ptr.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterContigs : public BaseFilter {

public:

  FilterContigs() = default;
  ~FilterContigs() override = default;

  [[nodiscard]] virtual std::unique_ptr<ContigDB> applyFilter(const ContigDB& contig) const = 0;
  [[nodiscard]] FilterBaseType filterType() const override { return FilterBaseType::CONTIG_FILTER; }

private:

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Base class for a filter at the OffsetDB level. Filters variants within an offset.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterOffsets : public BaseFilter {

public:

  FilterOffsets() = default;
  ~FilterOffsets() override = default;

  [[nodiscard]] virtual std::unique_ptr<OffsetDB> applyFilter(const OffsetDB& variant) const = 0;
  [[nodiscard]] FilterBaseType filterType() const override { return FilterBaseType::OFFSET_FILTER; }

private:

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Base class for a filter at the variant level. Filters individual variants.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FilterVariants : public BaseFilter {

public:

  FilterVariants() = default;
  ~FilterVariants() override = default;

  [[nodiscard]] virtual bool applyFilter(const Variant& variant) const = 0;
  [[nodiscard]] FilterBaseType filterType() const override { return FilterBaseType::VARIANT_FILTER; }

private:

};



} // Namespace

#endif //KGL_VARIANT_FILTER_TYPE_H
