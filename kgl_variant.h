///
// Created by kellerberrin on 13/10/17.
//

#ifndef KGL_VARIANT_H
#define KGL_VARIANT_H

#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_genome_types.h"
#include "kgl_alphabet_amino.h"
#include "kgl_genome_db.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class Variant; // Forward decl.
class ReadCountVariant; // Forward decl.
class SNPVariantDNA5; // Forward decl.
class CompoundVariant; // Forward decl.
class CompoundDelete; // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  virtual bool applyFilter(const Variant& variant) const = 0;
  virtual bool applyFilter(const ReadCountVariant& variant) const = 0;
  virtual bool applyFilter(const SNPVariantDNA5& variant) const = 0;
  virtual bool applyFilter(const CompoundDelete& variant) const = 0;

  virtual std::string filterName() const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class VariantSequenceType { UNKNOWN, CDS_CODING, INTRON, NON_CODING };

class VariantSequence {

public:

  VariantSequence(std::shared_ptr<const ContigFeatures> contig_ptr,
                ContigOffset_t contig_offset) : contig_ptr_(contig_ptr),
                                                contig_offset_(contig_offset),
                                                variant_genome_type_(VariantSequenceType::UNKNOWN) {}
  virtual ~VariantSequence() = default;

  std::string genomeOutput() const;  // Genome information text.

  std::shared_ptr<const ContigFeatures> contig() const { return contig_ptr_; }
  ContigOffset_t offset() const { return contig_offset_; }

  VariantSequenceType type() const { return genomeType(); };
  std::string typestr() const;
  const GeneVector& geneMembership() const { genomeType(); return gene_membership_; }

private:

  std::shared_ptr<const ContigFeatures> contig_ptr_;    // The contig.
  ContigOffset_t contig_offset_;                        // Location on the contig.
  VariantSequenceType variant_genome_type_;             // Non-coding, intron or coding.
  GeneVector gene_membership_;                          // Membership includes intron variants.

  VariantSequenceType genomeType(); // Lazy evaluation of gene membership and variant type.
  VariantSequenceType genomeType() const { return const_cast<VariantSequence*>(this)->genomeType(); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class Variant : public VariantSequence {

public:


  Variant(const std::shared_ptr<const ContigFeatures> contig_ptr,
          ContigOffset_t contig_offset) : VariantSequence(contig_ptr, contig_offset) {}
  Variant(const Variant& variant) = default;
  ~Variant() override = default;
  bool filterVariant(const VariantFilter& filter) const { return applyFilter(filter); }

  const ContigId_t& contigId() const { return contig()->contigId(); }
  ContigOffset_t contigOffset() const { return offset(); }
  bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };
  friend std::ostream& operator<<(std::ostream &os, const Variant& variant) { os << variant.output();  return os; }

  virtual std::string output() const = 0;
  virtual std::string mutation() const = 0;
  virtual bool isCompound() const { return false; }

private:

  virtual bool applyFilter(const VariantFilter& filter) const { return filter.applyFilter(*this); }
  virtual bool equivalent(const Variant& cmp_var) const = 0;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
