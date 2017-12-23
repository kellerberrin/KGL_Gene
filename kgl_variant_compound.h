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

  CompoundVariant(const std::string& variant_source,
                  std::shared_ptr<const ContigFeatures> contig_ptr,
                  ContigOffset_t contig_offset,
                  Phred_t quality,
                  const CompoundVariantMap& variant_map) : Variant(variant_source, contig_ptr, contig_offset, quality),
                                                           variant_map_(variant_map) {}
  ~CompoundVariant() override = default;

  size_t size() const override { return variant_map_.size(); }

  bool equivalent(const Variant& cmp_var) const override;

  const CompoundVariantMap& getMap() const { return variant_map_; }

  std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const;

protected:

  const CompoundVariantMap variant_map_;

  std::string location(char delimiter, VariantOutputIndex output_index) const;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound insert. Contiguous SNP insertions insert new nucleotides.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CompoundInsert : public CompoundVariant {

public:

  CompoundInsert(const std::string& variant_source,
                 std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 Phred_t quality,
                 const CompoundVariantMap& variant_map) : CompoundVariant(variant_source,
                                                                          contig_ptr,
                                                                          contig_offset,
                                                                          quality,
                                                                          variant_map) {}
  CompoundInsert(const CompoundInsert&) = default;
  ~CompoundInsert() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<CompoundInsert>(std::make_shared<CompoundInsert>(*this)); }

  VariantType variantType() const override { return VariantType::COMPOUND_INSERT; }

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound delete. Contiguous SNP deletions delete nucleotides.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CompoundDelete : public CompoundVariant {

public:

  CompoundDelete(const std::string& variant_source,
                 std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 Phred_t quality,
                 const CompoundVariantMap& variant_map) : CompoundVariant(variant_source,
                                                                          contig_ptr,
                                                                          contig_offset,
                                                                          quality,
                                                                          variant_map) {}
  CompoundDelete(const CompoundDelete&) = default;
  ~CompoundDelete() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<CompoundDelete>(std::make_shared<CompoundDelete>(*this)); }

  VariantType variantType() const override { return VariantType::COMPOUND_DELETE; }

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound SNP. where SNPs act on the same codon and therefore have a combined change to the AA.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CompoundSNP : public CompoundVariant {

public:

  CompoundSNP(const std::string& variant_source,
              std::shared_ptr<const ContigFeatures> contig_ptr,
              ContigOffset_t contig_offset,
              Phred_t quality,
              const CompoundVariantMap& variant_map) : CompoundVariant(variant_source,
                                                                       contig_ptr,
                                                                       contig_offset,
                                                                       quality,
                                                                       variant_map) {}
  CompoundSNP(const CompoundSNP&) = default;
  ~CompoundSNP() override = default;

  // Polymorphic copy constructor
  std::shared_ptr<Variant> clone() const override { return std::shared_ptr<CompoundSNP>(std::make_shared<CompoundSNP>(*this)); }

  VariantType variantType() const override { return VariantType::COMPOUND_SNP; }

  bool codonMutation( ContigOffset_t& codon_offset,
                      AminoAcid::Alphabet& reference_amino,
                      AminoAcid::Alphabet& mutant_amino) const;


private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                            ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                            SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_COMPOUND_H
