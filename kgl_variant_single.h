//
// Created by kellerberrin on 23/12/17.
//

#ifndef KGL_VARIANT_SINGLE_H
#define KGL_VARIANT_SINGLE_H


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

class SingleVariant : public Variant {

public:

  SingleVariant(const std::string& variant_source,
            std::shared_ptr<const ContigFeatures> contig_ptr,
            ContigOffset_t contig_offset,
            Phred_t quality,
            std::shared_ptr<const VariantEvidence> evidence_ptr,
            DNA5::Alphabet reference ) : Variant(variant_source, contig_ptr, contig_offset, quality),
                                           evidence_ptr_(evidence_ptr),
                                           reference_(reference) {}

  ~SingleVariant() override = default;

  size_t size() const override { return 1; }

  std::shared_ptr<const VariantEvidence> evidence() const { return evidence_ptr_; }

  DNA5::Alphabet reference() const { return reference_; }
  virtual char mutantChar() const = 0;

  std::string suboutput(char delimiter, VariantOutputIndex output_index) const;

  // complement base if -ve strand and coding or intron.
  CodingDNA5::Alphabet strandReference() const { return strandNucleotide(reference()); }


protected:

  CodingDNA5::Alphabet strandNucleotide(DNA5::Alphabet nucleotide) const;

private:

  const std::shared_ptr<const VariantEvidence> evidence_ptr_;
  DNA5::Alphabet reference_;

  std::string submutation(char delimiter, VariantOutputIndex output_index) const;

  std::string subname() const { return "S" + name(); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPVariant : public SingleVariant {

public:

  SNPVariant(const std::string& variant_source,
             std::shared_ptr<const ContigFeatures> contig_ptr,
             ContigOffset_t contig_offset,
             Phred_t quality,
             std::shared_ptr<const VariantEvidence> evidence_ptr,
             DNA5::Alphabet reference,
             DNA5::Alphabet mutant) : SingleVariant(variant_source,
                                                           contig_ptr,
                                                           contig_offset,
                                                           quality,
                                                           evidence_ptr,
                                                           reference), mutant_(mutant) {}

  SNPVariant(const SNPVariant& variant) = default;
  ~SNPVariant() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<SNPVariant>(std::make_shared<SNPVariant>(*this)); }

  VariantType variantType() const override { return VariantType::SNP; }

  bool equivalent(const Variant& cmp_var) const override;

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  CodingDNA5::Alphabet strandMutant() const { return strandNucleotide(mutant()); }

  DNA5::Alphabet mutant() const { return mutant_; }

  char mutantChar() const override { return DNA5::convertToChar(mutant()); }

  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

  bool codonMutation(ContigOffset_t &codon_offset,
                     ContigSize_t &base_in_codon,
                     AminoAcid::Alphabet &reference_amino,
                     AminoAcid::Alphabet &mutant_amino) const;

private:

  DNA5::Alphabet mutant_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple deletion variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DeleteVariant : public SingleVariant {

public:

  DeleteVariant(const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                ContigOffset_t contig_offset,
                Phred_t quality,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference) : SingleVariant(variant_source,
                                                              contig_ptr,
                                                              contig_offset,
                                                              quality,
                                                              evidence_ptr,
                                                              reference) {}

  DeleteVariant(const DeleteVariant& variant) = default;
  ~DeleteVariant() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<DeleteVariant>(std::make_shared<DeleteVariant>(*this)); }

  VariantType variantType() const override { return VariantType::DELETE; }

  bool equivalent(const Variant& cmp_var) const override;

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  virtual char mutantChar() const override { return '-'; }

  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple insertion variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InsertVariant : public SingleVariant {

public:

  InsertVariant(const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                ContigOffset_t contig_offset,
                Phred_t quality,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference,
                DNA5::Alphabet mutant) : SingleVariant(variant_source,
                                                       contig_ptr,
                                                       contig_offset,
                                                       quality,
                                                       evidence_ptr,
                                                       reference), mutant_(mutant) {}

  InsertVariant(const InsertVariant& variant) = default;
  ~InsertVariant() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<InsertVariant>(std::make_shared<InsertVariant>(*this)); }

  VariantType variantType() const override { return VariantType::INSERT; }

  bool equivalent(const Variant& cmp_var) const override;

  CodingDNA5::Alphabet strandMutant() const { return strandNucleotide(mutant()); }

  DNA5::Alphabet mutant() const { return mutant_; }

  virtual char mutantChar() const override { return DNA5::convertToChar(mutant()); }

  // This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

  DNA5::Alphabet mutant_;

};


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_SINGLE_H
