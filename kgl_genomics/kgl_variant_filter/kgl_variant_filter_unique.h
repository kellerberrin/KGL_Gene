//
// Created by kellerberrin on 17/07/23.
//

#ifndef KGL_VARIANT_FILTER_UNIQUE_H
#define KGL_VARIANT_FILTER_UNIQUE_H


#include "kgl_variant_filter_db_genome.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// For a haploid organism, such as the blood stage of P.Falciparum, only one variant can occur at a particular offset.
// These filters attempt to resolve the situation where there is more than 1 variant specified at a location.
// The filters assume that all variants have been converted to "canonical" format where the cigar format of an SNP is '1X',
// Deletes are '1MnD' and Inserts are '1MnI'.
//
// Note that an SNP and Indel specifying the same offset is not a problem, since  by convention, canonical
// indels actually occur at the next (+1) offset..
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter Unique variants for each offset.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Selection logic is somewhat involved because canonical INDEL variants actually operate on the NEXT (offset+1) offset.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This filter selects a variant from multiple variants at a particular offset randomly.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class RandomUniqueFilter : public FilterContigs {

public:

  RandomUniqueFilter() { filterName("RandomUniqueFilter"); }

  ~RandomUniqueFilter() override = default;

  [[nodiscard]] std::unique_ptr <ContigDB> applyFilter(const ContigDB &contig) const override { return filterUnique(contig); }
  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<RandomUniqueFilter>(*this); }


protected:

  // Selects a random unique variant by default.
  [[nodiscard]] virtual std::shared_ptr<const Variant>
  selectUnique(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const { return selectRandom(variant_vector); }

  [[nodiscard]] std::shared_ptr<const Variant> selectRandom(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  [[nodiscard]] std::shared_ptr<const Variant> selectFrequency(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const;

  [[nodiscard]] std::unique_ptr<ContigDB> filterUnique(const ContigDB &contig) const;
  [[nodiscard]] double getFrequency(const std::shared_ptr<const Variant>& variant_vector) const;
  void contigVector(std::unique_ptr<ContigDB>& contig_ptr, std::vector<std::shared_ptr<const Variant>>& offset_vector) const;

  constexpr static const char* AF_FIELD_{"AF"};

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This filter selects the most frequently occurring variant (highest probability of occurring) from multiple variants at an offset.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class FrequencyUniqueFilter : public RandomUniqueFilter {

public:

  FrequencyUniqueFilter() { filterName("FrequencyUniqueFilter"); }
  ~FrequencyUniqueFilter() override = default;

  [[nodiscard]] std::unique_ptr <ContigDB> applyFilter(const ContigDB &contig) const override { return filterUnique(contig); }
  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<FrequencyUniqueFilter>(*this); }


private:

  // Selects a random unique variant by default.
  [[nodiscard]] std::shared_ptr<const Variant>
  selectUnique(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const override { return selectFrequency(variant_vector); }


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This filter preferentially selects homozygous variants from multiple variants at an offset.
// If there is still ambiguity as to which variant modifies an offset then the most frequently
// occurring variant is selected.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class HomozygousUniqueFilter : public RandomUniqueFilter {

public:

  HomozygousUniqueFilter() { filterName("HomozygousUniqueFilter"); }
  ~HomozygousUniqueFilter() override = default;

  [[nodiscard]] std::unique_ptr <ContigDB> applyFilter(const ContigDB &contig) const override { return filterHomozygous(contig); }
  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<HomozygousUniqueFilter>(*this); }


private:

  // Preferentially selects homozygous variants.
  [[nodiscard]] std::unique_ptr <ContigDB> filterHomozygous(const ContigDB &contig) const;

};



} // Namespace.


#endif //KGL_VARIANT_FILTER_UNIQUE_H
