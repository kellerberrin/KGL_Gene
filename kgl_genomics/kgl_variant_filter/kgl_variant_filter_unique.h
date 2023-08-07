//
// Created by kellerberrin on 17/07/23.
//

#ifndef KGL_VARIANT_FILTER_UNIQUE_H
#define KGL_VARIANT_FILTER_UNIQUE_H


#include "kgl_variant_filter_db.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// For a haploid organism, such as the blood stage of P.Falciparum, only one variant can occur at a particular offset.
// These filters attempt to resolve the situation where there is more than 1 valid variant specified at a location.
// The filters assume that all variants have been converted to "canonical" format where the cigar format of an SNP is '1X',
// Deletes are '1MnD' and Inserts are '1MnI'.
//
// Note that an SNP and Indel specifying the same offset is not a problem, since  by convention, canonical
// Indels actually occur at the next (+1) offset..
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Unique variants for each offset.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Selection logic is somewhat involved because canonical INDEL variants actually operate on the NEXT (offset+1) offset.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RandomUniqueFilter : public FilterContigs {

public:

  RandomUniqueFilter() { filterName("RandomUniqueFilter"); }

  ~RandomUniqueFilter() override = default;

  [[nodiscard]] std::unique_ptr <ContigDB> applyFilter(const ContigDB &contig) const override { return filterUnique(contig); }
  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<RandomUniqueFilter>(*this); }


protected:

  // Selects a random unique variant by default.
  [[nodiscard]] virtual std::shared_ptr<const Variant> selectUnique(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const;
  [[nodiscard]] std::unique_ptr<ContigDB> filterUnique(const ContigDB &contig) const;
  void contigVector(std::unique_ptr<ContigDB>& contig_ptr, std::vector<std::shared_ptr<const Variant>>& offset_vector) const;

};


class FrequencyUniqueFilter : public RandomUniqueFilter {

public:

  FrequencyUniqueFilter() { filterName("FrequencyUniqueFilter"); }
  ~FrequencyUniqueFilter() override = default;

  [[nodiscard]] std::unique_ptr <ContigDB> applyFilter(const ContigDB &contig) const override { return filterUnique(contig); }
  [[nodiscard]] std::shared_ptr <BaseFilter> clone() const override { return std::make_shared<FrequencyUniqueFilter>(*this); }


private:

  // Selects a random unique variant by default.
  [[nodiscard]] std::shared_ptr<const Variant> selectUnique(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const override;
  [[nodiscard]] double getFrequency(const std::shared_ptr<const Variant>& variant_vector) const;


  constexpr static const char* AF_FIELD_{"AF"};

};



} // Namespace.


#endif //KGL_VARIANT_FILTER_UNIQUE_H
