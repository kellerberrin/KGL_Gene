//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_COMPOUND_H
#define KGL_VARIANT_COMPOUND_H


#include "kgl_variant.h"
#include "kgl_variant_snp.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound variant. A collection of feature aligned and contiguous variants. Insertions and Deletions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using CompoundVariantMap = std::map<ContigOffset_t, std::shared_ptr<const SubordinateSNP>>;

class CompoundVariant : public Variant {

public:

  CompoundVariant(const std::string& variant_source,
                  std::shared_ptr<const ContigFeatures> contig_ptr,
                  ContigOffset_t contig_offset,
                  const CompoundVariantMap& variant_map) : Variant(variant_source, contig_ptr, contig_offset),
                                                           variant_map_(variant_map) {}
  ~CompoundVariant() override = default;

  bool isCompound() const override { return true; }

  bool equivalent(const Variant& cmp_var) const override;

  const CompoundVariantMap& getMap() const { return variant_map_; }

protected:

  const CompoundVariantMap variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound insert. Contiguous SNP insertions insert new nucleotides.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CompoundInsert : public CompoundVariant {

public:

  CompoundInsert(const std::string& variant_source,
                 std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 const CompoundVariantMap& variant_map) : CompoundVariant(variant_source,
                                                                          contig_ptr,
                                                                          contig_offset,
                                                                          variant_map) {}
  ~CompoundInsert() override = default;

  static constexpr const char* VARIANT_NAME = "CINS";
  std::string name() const override { return VARIANT_NAME; }

private:


  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string output(char delimter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
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
                 const CompoundVariantMap& variant_map) : CompoundVariant(variant_source,
                                                                          contig_ptr,
                                                                          contig_offset,
                                                                          variant_map) {}
  ~CompoundDelete() override = default;

  static constexpr const char* VARIANT_NAME = "CDEL";
  std::string name() const override { return VARIANT_NAME; }

private:

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string output(char delimter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
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
              const CompoundVariantMap& variant_map) : CompoundVariant(variant_source,
                                                                       contig_ptr,
                                                                       contig_offset,
                                                                       variant_map) {}
  ~CompoundSNP() override = default;

  static constexpr const char* VARIANT_NAME = "CSNP";
  std::string name() const override { return VARIANT_NAME; }

private:


  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string output(char delimiter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;
  bool codonMutation( ContigOffset_t& codon_offset,
                      AminoAcid::Alphabet& reference_amino,
                      AminoAcid::Alphabet& mutant_amino) const;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_COMPOUND_H
