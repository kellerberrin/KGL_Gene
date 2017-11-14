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
                 Nucleotide_DNA5_t reference,
                 Nucleotide_DNA5_t mutant)
  : ReadCountVariant(contig_ptr, contig_offset, read_count, mutant_count, count_array, count_array_size),
    reference_(reference),
    mutant_(mutant) {}

  SNPVariantDNA5(const SNPVariantDNA5& variant) = default;
  ~SNPVariantDNA5() override = default;

  bool equivalent(const Variant& cmp_var) const override;

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;

  Nucleotide_DNA5_t reference() const { return reference_; }
  Nucleotide_ExtendedDNA5 mutant() const { return mutant_; }

  // complement base if -ve strand and coding or intron.
  Nucleotide_DNA5_t strandReference() const { return strandNucleotide(reference()); }
  // complement base if -ve strand and coding or intron.
  Nucleotide_ExtendedDNA5 strandMutant() const { return strandNucleotide(mutant()); }

  std::string output(char delimiter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

private:

  Nucleotide_DNA5_t reference_;
  Nucleotide_ExtendedDNA5 mutant_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  Nucleotide_ExtendedDNA5 strandNucleotide(Nucleotide_ExtendedDNA5 nucleotide) const;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_SNP_H
