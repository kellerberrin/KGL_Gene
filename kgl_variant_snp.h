//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_SNP_H
#define KGL_VARIANT_SNP_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_genome_types.h"
#include "kgl_alphabet_amino.h"
#include "kgl_variant.h"
#include "kgl_genome_db.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A read count variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ReadCountVariant : public Variant {

public:

  ReadCountVariant(const std::shared_ptr<const ContigFeatures> contig_ptr,
                   ContigOffset_t contig_offset,
                   NucleotideReadCount_t read_count,
                   NucleotideReadCount_t mutant_count,
                   NucleotideReadCount_t const count_array[],
                   ContigSize_t count_array_size) : Variant(contig_ptr, contig_offset),
                                                    read_count_(read_count),
                                                    mutant_count_(mutant_count) {

    for(ContigOffset_t idx = 0; idx < count_array_size; ++idx) {

      count_array_.push_back(count_array[idx]);

    }

  }
  ReadCountVariant(const ReadCountVariant& variant) = default;
  ~ReadCountVariant() override = default;

  NucleotideReadCount_t readCount() const { return read_count_; }
  NucleotideReadCount_t mutantCount() const { return mutant_count_; }
  const std::vector<NucleotideReadCount_t>& countArray() const { return count_array_; }

  double proportion() const { return static_cast<double>(mutant_count_) / static_cast<double>(read_count_); }

private:

  NucleotideReadCount_t read_count_;
  NucleotideReadCount_t mutant_count_;
  std::vector<NucleotideReadCount_t> count_array_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPVariantDNA5 : public ReadCountVariant {

public:

  SNPVariantDNA5(const std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 NucleotideReadCount_t read_count,
                 NucleotideReadCount_t mutant_count,
                 NucleotideReadCount_t const count_array[],
                 ContigSize_t  count_array_size,
                 typename NucleotideColumn_DNA5::NucleotideType reference,
                 typename NucleotideColumn_DNA5::NucleotideType mutant)
  : ReadCountVariant(contig_ptr, contig_offset, read_count, mutant_count, count_array, count_array_size),
    reference_(reference),
    mutant_(mutant) {}

  SNPVariantDNA5(const SNPVariantDNA5& variant) = default;

  ~SNPVariantDNA5() override = default;

  bool equivalent(const Variant& cmp_var) const override;

  const typename NucleotideColumn_DNA5::NucleotideType& reference() const { return reference_; }
  const typename NucleotideColumn_DNA5::NucleotideType& mutant() const { return mutant_; }

  std::string output(char delimiter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

private:

  typename NucleotideColumn_DNA5::NucleotideType reference_;
  typename NucleotideColumn_DNA5::NucleotideType mutant_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_SNP_H
