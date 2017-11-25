//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_VARIANT_COMPOUND_H
#define KGL_VARIANT_COMPOUND_H


#include "kgl_variant.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound variant. A collection of feature aligned and contiguous variants. Insertions and Deletions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CompoundVariant : public Variant {

public:

  CompoundVariant(const std::string& variant_source,
                  std::shared_ptr<const ContigFeatures> contig_ptr,
                  ContigOffset_t contig_offset) : Variant(variant_source, contig_ptr, contig_offset) {}
  ~CompoundVariant() override = default;

  bool isCompound() const override { return true; }

  bool equivalent(const Variant& cmp_var) const override;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound insert. Contiguous SNP insertions insert new nucleotides.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using CompoundVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;

class CompoundInsert : public CompoundVariant {

public:

  CompoundInsert(const std::string& variant_source,
                 std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 const CompoundVariantMap& variant_map) : CompoundVariant(variant_source, contig_ptr, contig_offset),
                                                         variant_map_(variant_map) {}
  ~CompoundInsert() override = default;

  const CompoundVariantMap& getMap() const { return variant_map_; }

  static constexpr const char* VARIANT_TYPE = "COMPOUND_INSERT";

private:

  const CompoundVariantMap variant_map_;

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
                 const CompoundVariantMap& variant_map) : CompoundVariant(variant_source, contig_ptr, contig_offset),
                                                         variant_map_(variant_map) {}
  ~CompoundDelete() override = default;

  const CompoundVariantMap& getMap() const { return variant_map_; }

  static constexpr const char* VARIANT_TYPE = "COMPOUND_DELETE";

private:

  const CompoundVariantMap variant_map_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string output(char delimter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound SNP. where SNPs act on the same codon and therefore have a combined change to the AA.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using SNPCompoundVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;

class CompoundSNP : public CompoundVariant {

public:

  CompoundSNP(const std::string& variant_source,
              std::shared_ptr<const ContigFeatures> contig_ptr,
              ContigOffset_t contig_offset,
              const SNPCompoundVariantMap& variant_map) : CompoundVariant(variant_source, contig_ptr, contig_offset),
                                                         variant_map_(variant_map) {}
  ~CompoundSNP() override = default;

  const SNPCompoundVariantMap& getMap() const { return variant_map_; }

  static constexpr const char* VARIANT_TYPE = "COMPOUND_SNP";

private:

  const SNPCompoundVariantMap variant_map_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string output(char delimiter, VariantOutputIndex output_index) const override;
  std::string mutation(char delimiter, VariantOutputIndex output_index) const override;
  bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const override;


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_COMPOUND_H
