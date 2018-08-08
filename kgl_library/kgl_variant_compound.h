//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_COMPOUND_H
#define KGL_VARIANT_COMPOUND_H


#include "kgl_variant.h"
#include "kgl_variant_single.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound variant. A collection of feature aligned and contiguous variants. Insertions and Deletions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using CompoundVariantMap = std::map<ContigOffset_t, std::shared_ptr<const SingleVariant>>;

class CompoundVariant : public Variant {

public:

  CompoundVariant(const GenomeId_t& genome_id,
                  const ContigId_t& contig_id,
                  PhaseId_t phase_id,
                  ContigOffset_t contig_offset,
                  std::shared_ptr<const VariantEvidence> evidence_ptr,
                  const CompoundVariantMap& variant_map) : Variant(genome_id, contig_id, phase_id, contig_offset, evidence_ptr),
                                                           variant_map_(variant_map) {}
  ~CompoundVariant() override = default;

  size_t size() const override { return variant_map_.size(); }

  void updatePhaseId(PhaseId_t phase_id) override;

  bool equivalent(const Variant& cmp_var) const override;
  bool lessThan(const Variant& cmp_var) const override;

  const CompoundVariantMap& getMap() const { return variant_map_; }

  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const;

protected:

  const CompoundVariantMap variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound insert. Contiguous SNP insertions insert new nucleotides.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CompoundInsert : public CompoundVariant {

public:

  CompoundInsert(const GenomeId_t& genome_id,
                 const ContigId_t& contig_id,
                 PhaseId_t phase_id,
                 ContigOffset_t contig_offset,
                 std::shared_ptr<const VariantEvidence> evidence_ptr,
                 const CompoundVariantMap& variant_map) : CompoundVariant(genome_id,
                                                                          contig_id,
                                                                          phase_id,
                                                                          contig_offset,
                                                                          evidence_ptr,
                                                                          variant_map) {}
  CompoundInsert(const CompoundInsert&) = default;
  ~CompoundInsert() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<CompoundInsert>(std::make_shared<CompoundInsert>(*this)); }

  VariantType variantType() const override { return VariantType::COMPOUND_INSERT; }

  bool mutateSequence(SignedOffset_t offset_adjust,
                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                      SignedOffset_t& sequence_size_modify) const override;

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound delete. Contiguous SNP deletions delete nucleotides.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CompoundDelete : public CompoundVariant {

public:

  CompoundDelete(const GenomeId_t& genome_id,
                 const ContigId_t& contig_id,
                 PhaseId_t phase_id,
                 ContigOffset_t contig_offset,
                 std::shared_ptr<const VariantEvidence> evidence_ptr,
                 const CompoundVariantMap& variant_map) : CompoundVariant(genome_id,
                                                                          contig_id,
                                                                          phase_id,
                                                                          contig_offset,
                                                                          evidence_ptr,
                                                                          variant_map) {}
  CompoundDelete(const CompoundDelete&) = default;
  ~CompoundDelete() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<CompoundDelete>(std::make_shared<CompoundDelete>(*this)); }

  VariantType variantType() const override { return VariantType::COMPOUND_DELETE; }

  bool mutateSequence(SignedOffset_t offset_adjust,
                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                      SignedOffset_t& sequence_size_modify) const override;


private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_COMPOUND_H
