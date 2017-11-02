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

using CompoundVariantMap = std::map<ContigOffset_t, std::shared_ptr<const Variant>>;
class CompoundVariant : public Variant {

public:

  CompoundVariant(std::shared_ptr<const ContigFeatures> contig_ptr,
                  ContigOffset_t contig_offset,
                  const CompoundVariantMap& variant_map) : Variant(contig_ptr, contig_offset),
                                                           variant_map_(variant_map) {}
  ~CompoundVariant() override = default;

  bool equivalent(const Variant& cmp_var) const override;

protected:

  const CompoundVariantMap variant_map_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound delete. 3 contiguous SNP deletions aligned on a Gene codon boundary delete a single Amino Acid
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CompoundDelete : public CompoundVariant {

public:

  CompoundDelete(std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 const CompoundVariantMap variant_map) : CompoundVariant(contig_ptr, contig_offset, variant_map) {}
  ~CompoundDelete() override = default;

private:


  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }
  std::string output() const;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_COMPOUND_H
