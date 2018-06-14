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

  SingleVariant(const GenomeId_t& genome_id,
                const ContigId_t& contig_id,
                PhaseId_t phase_id,
                ContigOffset_t contig_offset,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference ) : Variant(genome_id, contig_id, phase_id, contig_offset, evidence_ptr),
                                             reference_(reference) {}

  ~SingleVariant() override = default;

  size_t size() const override { return 1; }

  DNA5::Alphabet reference() const { return reference_; }

  virtual char mutantChar() const = 0;

  std::string suboutput(char delimiter, VariantOutputIndex output_index, bool detail) const;

private:

  DNA5::Alphabet reference_;

  std::string submutation(char delimiter, VariantOutputIndex output_index) const;

  std::string subname() const { return "S" + name(); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPVariant : public SingleVariant {

public:

  SNPVariant(const GenomeId_t& genome_id,
             const ContigId_t& contig_id,
             PhaseId_t phase_id,
             ContigOffset_t contig_offset,
             std::shared_ptr<const VariantEvidence> evidence_ptr,
             DNA5::Alphabet reference,
             DNA5::Alphabet mutant) : SingleVariant(genome_id,
                                                    contig_id,
                                                    phase_id,
                                                    contig_offset,
                                                    evidence_ptr,
                                                    reference), mutant_(mutant) {}

  SNPVariant(const SNPVariant& variant) = default;
  ~SNPVariant() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<SNPVariant>(std::make_shared<SNPVariant>(*this)); }

  VariantType variantType() const override { return VariantType::SNP; }

  bool equivalent(const Variant& cmp_var) const override;
  bool lessThan(const Variant& cmp_var) const override;

  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  bool mutateSequence(SignedOffset_t offset_adjust,
                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                      SignedOffset_t& sequence_size_modify) const override;

  DNA5::Alphabet mutant() const { return mutant_; }

  char mutantChar() const override { return DNA5::convertToChar(mutant()); }

  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

private:

  DNA5::Alphabet mutant_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple deletion variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DeleteVariant : public SingleVariant {

public:

  DeleteVariant(const GenomeId_t& genome_id,
                const ContigId_t& contig_id,
                PhaseId_t phase_id,
                ContigOffset_t contig_offset,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference) : SingleVariant(genome_id,
                                                          contig_id,
                                                          phase_id,
                                                          contig_offset,
                                                          evidence_ptr,
                                                          reference) {}

  DeleteVariant(const DeleteVariant& variant) = default;
  ~DeleteVariant() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<DeleteVariant>(std::make_shared<DeleteVariant>(*this)); }

  VariantType variantType() const override { return VariantType::DELETE; }

  bool equivalent(const Variant& cmp_var) const override;
  bool lessThan(const Variant& cmp_var) const override;

  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  bool mutateSequence(SignedOffset_t offset_adjust,
                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                      SignedOffset_t& sequence_size_modify) const override;

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

  InsertVariant(const GenomeId_t& genome_id,
                const ContigId_t& contig_id,
                PhaseId_t phase_id,
                ContigOffset_t contig_offset,
                std::shared_ptr<const VariantEvidence> evidence_ptr,
                DNA5::Alphabet reference,
                DNA5::Alphabet mutant) : SingleVariant(genome_id,
                                                       contig_id,
                                                       phase_id,
                                                       contig_offset,
                                                       evidence_ptr,
                                                       reference), mutant_(mutant) {}

  InsertVariant(const InsertVariant& variant) = default;
  ~InsertVariant() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<InsertVariant>(std::make_shared<InsertVariant>(*this)); }

  VariantType variantType() const override { return VariantType::INSERT; }

  bool equivalent(const Variant& cmp_var) const override;
  bool lessThan(const Variant& cmp_var) const override;

  DNA5::Alphabet mutant() const { return mutant_; }

  virtual char mutantChar() const override { return DNA5::convertToChar(mutant()); }

  bool mutateSequence(SignedOffset_t offset_adjust,
                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                      SignedOffset_t& sequence_size_modify) const override;


  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

  DNA5::Alphabet mutant_;

};


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_SINGLE_H
