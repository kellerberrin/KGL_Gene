///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant_single.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to coding sequences only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CodingFilter : public VariantFilter {

public:

  explicit CodingFilter() {}
  ~CodingFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  bool implementFilter(const Variant& variant) const;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants to SNPs (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPFilter : public VariantFilter {

public:

  explicit SNPFilter() {}
  ~SNPFilter() override = default;

  std::string filterName() const final { return "Filter Single SNP and Compound MNP (exclude indels)"; }

  bool applyFilter(const SNPVariant&) const override { return true; }
  bool applyFilter(const DeleteVariant&) const override { return false; }
  bool applyFilter(const InsertVariant&) const override { return false; }
  bool applyFilter(const CompoundDelete&) const override { return false; }
  bool applyFilter(const CompoundInsert&) const override { return false; }
  bool applyFilter(const CompoundSNP&) const override { return true; }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter to Delete variants (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DeleteFilter : public VariantFilter {

public:

  explicit DeleteFilter() {}
  ~DeleteFilter() override = default;

  std::string filterName() const final { return "Filter Single and Compound Delete Variants"; }

  bool applyFilter(const SNPVariant&) const override { return false; }
  bool applyFilter(const DeleteVariant&) const override { return true; }
  bool applyFilter(const InsertVariant&) const override { return false; }
  bool applyFilter(const CompoundDelete&) const override { return true; }
  bool applyFilter(const CompoundInsert&) const override { return false; }
  bool applyFilter(const CompoundSNP&) const override { return false; }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter to Insert variants (single and compound)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InsertFilter : public VariantFilter {

public:

  explicit InsertFilter() {}
  ~InsertFilter() override = default;

  std::string filterName() const final { return "Filter Single and Compound Insert Variants"; }

  bool applyFilter(const SNPVariant&) const override { return false; }
  bool applyFilter(const DeleteVariant&) const override { return false; }
  bool applyFilter(const InsertVariant&) const override { return true; }
  bool applyFilter(const CompoundDelete&) const override { return false; }
  bool applyFilter(const CompoundInsert&) const override { return true; }
  bool applyFilter(const CompoundSNP&) const override { return false; }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter variants by quality (-10log10 {prob variant is in error})
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class QualityFilter : public VariantFilter {

public:

  explicit QualityFilter(Phred_t quality) : quality_(quality) {}
  ~QualityFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const DeleteVariant& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const InsertVariant& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const CompoundDelete& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const CompoundInsert& variant) const override { return variant.quality() >= quality_; }
  bool applyFilter(const CompoundSNP& variant) const override { return variant.quality() >= quality_; }

private:

  const Phred_t quality_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {}
  ~ContigFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  const ContigId_t contig_ident_;

  bool implementFilter(const Variant& variant) const;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular gene.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneFilter : public VariantFilter {

public:

  explicit GeneFilter(const FeatureIdent_t& gene_ident) : gene_ident_(gene_ident) {}
  ~GeneFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  const FeatureIdent_t gene_ident_;

  bool implementFilter(const Variant& variant) const;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular CDS sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SequenceFilter : public VariantFilter {

public:

  explicit SequenceFilter(const FeatureIdent_t& sequence_ident) : sequence_ident_(sequence_ident) {}
  ~SequenceFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  const FeatureIdent_t sequence_ident_;

  bool implementFilter(const Variant& variant) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remove synonymous coding SNPs and compound snps.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SynonymousFilter : public VariantFilter {

public:

  explicit SynonymousFilter() = default;
  ~SynonymousFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  bool implementFilter(const Variant& variant) const;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T, typename... Args>
class NotFilter : public VariantFilter {

public:

  explicit NotFilter(Args... args) : filter_ptr_(std::make_unique<T>(args...)) {}
  ~NotFilter() override = default;

  bool applyFilter(const SNPVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const DeleteVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const InsertVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return not filter_ptr_->applyFilter(variant); }


  std::string filterName() const final { return "NOT (" + filter_ptr_->filterName() + ")"; }

private:

  std::unique_ptr<const T> filter_ptr_;

};

using NotCodingFilter = NotFilter<CodingFilter>;
using NotContigFilter = NotFilter<ContigFilter>;
using NotGeneFilter = NotFilter<GeneFilter>;
using NotSequenceFilter = NotFilter<SequenceFilter>;
using NotSynonymousFilter = NotFilter<SynonymousFilter>;
using NotSNPFilter = NotFilter<SNPFilter>;
using NotDeleteFilter = NotFilter<DeleteFilter>;
using NotInsertFilter = NotFilter<InsertFilter>;



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
