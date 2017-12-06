///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant_snp.h"
#include "kgl_variant_compound.h"
#include "kgl_genome_db.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the minimum read count SNP generation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ReadCountFilter : public VariantFilter {

public:

  explicit ReadCountFilter(NucleotideReadCount_t read_count) : read_count_(read_count) {}
  ~ReadCountFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return true; }  // default
  bool applyFilter(const CompoundInsert& variant) const override { return true; }
  bool applyFilter(const CompoundSNP& variant) const override { return true; }

  std::string filterName() const final;

private:

  bool implementFilter(const SNPVariant& variant) const;

  NucleotideReadCount_t read_count_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the minimum mutant read proportion in a candidate SNP.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MutantProportionFilter : public VariantFilter {

public:

  explicit MutantProportionFilter(double proportion) : mutant_proportion_(proportion) {}
  ~MutantProportionFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return true; }  // default
  bool applyFilter(const CompoundInsert& variant) const override { return true; }
  bool applyFilter(const CompoundSNP& variant) const override { return true; }

  std::string filterName() const final;

private:

  bool implementFilter(const SNPVariant& variant) const;

  double mutant_proportion_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to coding sequences only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CodingFilter : public VariantFilter {

public:

  explicit CodingFilter() {}
  ~CodingFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  bool implementFilter(const Variant& variant) const;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident) : contig_ident_(contig_ident) {}
  ~ContigFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
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

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
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

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
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

class SynonymousSNPFilter : public VariantFilter {

public:

  explicit SynonymousSNPFilter() = default;
  ~SynonymousSNPFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const SNPVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return implementFilter(variant); }

private:

  bool implementFilter(const Variant& variant) const;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter deletion SNPs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DeleteSNPFilter : public VariantFilter {

public:

  DeleteSNPFilter() = default;
  ~DeleteSNPFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariant& variant) const override {

    return variant.variantType() == VariantType::DELETE;

  }
  bool applyFilter(const CompoundDelete& variant) const override { return true; }
  bool applyFilter(const CompoundInsert& variant) const override { return false; }
  bool applyFilter(const CompoundSNP& variant) const override { return false; }


  std::string filterName() const final { return "Delete nucleotide base '-' SNPs only"; }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter mutant SNPs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MutantSNPFilter : public VariantFilter {

public:

  MutantSNPFilter() = default;
  ~MutantSNPFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariant& variant) const override {

    return variant.variantType() == VariantType::SNP;

  }
  bool applyFilter(const CompoundDelete& variant) const override { return false; }
  bool applyFilter(const CompoundInsert& variant) const override { return false; }
  bool applyFilter(const CompoundSNP& variant) const override { return true; }

  std::string filterName() const final { return "Mutation SNPs (no insert '+' or delete '-' SNPS)"; }

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Negation Filter
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T, typename... Args>
class NotFilter : public VariantFilter {

public:

  explicit NotFilter(Args... args) : filter_ptr_(std::make_unique<T>(args...)) {}
  ~NotFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const SNPVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundDelete& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundInsert& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CompoundSNP& variant) const override { return not filter_ptr_->applyFilter(variant); }


  std::string filterName() const final { return "NOT (" + filter_ptr_->filterName() + ")"; }

private:

  std::unique_ptr<const T> filter_ptr_;

};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
