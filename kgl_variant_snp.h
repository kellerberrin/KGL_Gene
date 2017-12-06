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
#include "kgl_variant_evidence.h"
#include "kgl_variant.h"
#include "kgl_genome_db.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A virtual class held in compound variants, produces a modified text output for coding variants.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SubordinateSNP : public Variant {

public:

  SubordinateSNP(const std::string& variant_source,
                 std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 std::shared_ptr<const VariantEvidence> evidence_ptr,
                 DNA5::Alphabet reference,
                 ExtendDNA5::Alphabet mutant) : Variant(variant_source, contig_ptr, contig_offset),
                                                evidence_ptr_(evidence_ptr),
                                                reference_(reference),
                                                mutant_(mutant) {}

  size_t size() const override { return 1; }

  VariantType variantType() const override;

  std::shared_ptr<const VariantEvidence> evidence() const { return evidence_ptr_; }

  DNA5::Alphabet reference() const { return reference_; }
  ExtendDNA5::Alphabet mutant() const { return mutant_; }

  std::string suboutput(char delimiter, VariantOutputIndex output_index) const;

private:

  const std::shared_ptr<const VariantEvidence> evidence_ptr_;
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
             std::shared_ptr<const ContigFeatures> contig_ptr,
             ContigOffset_t contig_offset,
             std::shared_ptr<const VariantEvidence> evidence_ptr,
             DNA5::Alphabet reference,
             ExtendDNA5::Alphabet mutant) : SubordinateSNP(variant_source,
                                                           contig_ptr,
                                                           contig_offset,
                                                           evidence_ptr,
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

  bool codonMutation(ContigOffset_t &codon_offset,
                     ContigSize_t &base_in_codon,
                     AminoAcid::Alphabet &reference_amino,
                     AminoAcid::Alphabet &mutant_amino) const;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  CodingDNA5::Alphabet strandNucleotide(DNA5::Alphabet nucleotide) const;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_SNP_H
