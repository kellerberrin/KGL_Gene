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

  ReadCountVariant(const std::string& variant_source,
                   const std::shared_ptr<const ContigFeatures> contig_ptr,
                   ContigOffset_t contig_offset,
                   NucleotideReadCount_t read_count,
                   NucleotideReadCount_t mutant_count,
                   NucleotideReadCount_t const count_array[],
                   ContigSize_t count_array_size) : Variant(variant_source, contig_ptr, contig_offset),
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
//  A virtual class held in compound variants, produces a modified text output for coding variants.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SubordinateSNP : public ReadCountVariant {

public:

  SubordinateSNP(const std::string& variant_source,
                const std::shared_ptr<const ContigFeatures> contig_ptr,
                ContigOffset_t contig_offset,
                NucleotideReadCount_t read_count,
                NucleotideReadCount_t mutant_count,
                NucleotideReadCount_t const count_array[],
                ContigSize_t  count_array_size,
                DNA5::Alphabet reference,
                ExtendDNA5::Alphabet mutant) : ReadCountVariant(variant_source,
                                                                contig_ptr,
                                                                contig_offset,
                                                                read_count,
                                                                mutant_count,
                                                                count_array,
                                                                count_array_size),
                                              reference_(reference),
                                              mutant_(mutant) {}

  size_t size() const override { return 1; }

  VariantType variantType() const override;

  DNA5::Alphabet reference() const { return reference_; }
  ExtendDNA5::Alphabet mutant() const { return mutant_; }


  std::string suboutput(char delimiter, VariantOutputIndex output_index) const;

private:

  DNA5::Alphabet reference_;
  ExtendDNA5::Alphabet mutant_;

  std::string submutation(char delimiter, VariantOutputIndex output_index) const;

  std::string subname() const { return "S" + name(); }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPVariant : public SubordinateSNP {

public:

  SNPVariant(const std::string& variant_source,
             const std::shared_ptr<const ContigFeatures> contig_ptr,
             ContigOffset_t contig_offset,
             NucleotideReadCount_t read_count,
             NucleotideReadCount_t mutant_count,
             NucleotideReadCount_t const count_array[],
             ContigSize_t  count_array_size,
             DNA5::Alphabet reference,
             ExtendDNA5::Alphabet mutant) : SubordinateSNP(variant_source,
                                                           contig_ptr,
                                                           contig_offset,
                                                           read_count,
                                                           mutant_count,
                                                           count_array,
                                                           count_array_size,
                                                           reference,
                                                           mutant) {}

  SNPVariant(const SNPVariant& variant) = default;
  ~SNPVariant() override = default;

  bool equivalent(const Variant& cmp_var) const override;

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;

  // complement base if -ve strand and coding or intron.
  CodingDNA5::Alphabet strandReference() const { return strandNucleotide(reference()); }
  CodingDNA5::Alphabet strandMutant() const { return strandNucleotide(ExtendDNA5::extendToBase(mutant())); }

  std::string output(char delimiter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

  bool SNPMutation(ContigOffset_t& codon_offset,
                   ContigSize_t& base_in_codon,
                   AminoAcid::Alphabet& reference_amino,
                   AminoAcid::Alphabet& mutant_amino) const;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  CodingDNA5::Alphabet strandNucleotide(DNA5::Alphabet nucleotide) const;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_SNP_H
