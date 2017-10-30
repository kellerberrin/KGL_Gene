///
// Created by kellerberrin on 16/10/17.
//

#ifndef KGL_FILTER_H
#define KGL_FILTER_H

#include "kgl_variant.h"
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
  bool applyFilter(const ReadCountVariant& variant) const override { return variant.readCount() >= read_count_; }
  bool applyFilter(const SNPVariantDNA5& variant) const override { return variant.readCount() >= read_count_; }
  bool applyFilter(const CodonDelete& variant) const override { return true; }  // default

  std::string filterName() const final;

private:

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
  bool applyFilter(const ReadCountVariant& variant) const override { return variant.proportion() >= mutant_proportion_; }
  bool applyFilter(const SNPVariantDNA5& variant) const override { return variant.proportion() >= mutant_proportion_; }
  bool applyFilter(const CodonDelete& variant) const override { return true; }  // default

  std::string filterName() const final;

private:

  double mutant_proportion_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to coding sequences only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InCDSFilter : public VariantFilter {

public:

  explicit InCDSFilter(std::shared_ptr<const GenomeDatabase> genome_db_ptr) : genome_db_ptr_(genome_db_ptr) {}
  ~InCDSFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const ReadCountVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const SNPVariantDNA5& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CodonDelete& variant) const override { return implementFilter(variant); }

private:

  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

  bool implementFilter(const Variant& variant) const;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular contig.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public VariantFilter {

public:

  explicit ContigFilter(const ContigId_t& contig_ident,
                      const std::shared_ptr<const GenomeDatabase> genome_db_ptr) : contig_ident_(contig_ident),
                                                                                   genome_db_ptr_(genome_db_ptr) {}
  ~ContigFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const ReadCountVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const SNPVariantDNA5& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CodonDelete& variant) const override { return implementFilter(variant); }

private:

  const ContigId_t contig_ident_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

  bool implementFilter(const Variant& variant) const;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter SNPs to a particular gene.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneFilter : public VariantFilter {

public:

  explicit GeneFilter(const FeatureIdent_t& gene_ident,
                      const std::shared_ptr<const GenomeDatabase> genome_db_ptr) : gene_ident_(gene_ident),
                                                                                   genome_db_ptr_(genome_db_ptr) {}
  ~GeneFilter() override = default;

  std::string filterName() const final;

  bool applyFilter(const Variant& variant) const override { return implementFilter(variant); } // redirect
  bool applyFilter(const ReadCountVariant& variant) const override { return implementFilter(variant); }
  bool applyFilter(const SNPVariantDNA5& variant) const override { return implementFilter(variant); }
  bool applyFilter(const CodonDelete& variant) const override { return implementFilter(variant); }

private:

  const FeatureIdent_t gene_ident_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

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
  bool applyFilter(const ReadCountVariant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariantDNA5& variant) const override {

    return variant.mutant() == NucleotideColumn_DNA5::DELETE_NUCLEOTIDE;

  }
  bool applyFilter(const CodonDelete& variant) const override { return true; }


  std::string filterName() const final { return "Delete nucleotide base '-' SNPs only"; }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter insert SNPs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class InsertSNPFilter : public VariantFilter {

public:

  InsertSNPFilter() = default;
  ~InsertSNPFilter() override = default;

  bool applyFilter(const Variant& variant) const override { return true; }  // default
  bool applyFilter(const ReadCountVariant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariantDNA5& variant) const override {

    return variant.mutant() == NucleotideColumn_DNA5::INSERT_SEQUENCE;

  }
  bool applyFilter(const CodonDelete& variant) const override { return true; }

  std::string filterName() const final { return "Insert nucleotide bases '+' SNPs only"; }

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
  bool applyFilter(const ReadCountVariant& variant) const override { return true; }  // default
  bool applyFilter(const SNPVariantDNA5& variant) const override {

    return variant.mutant() != NucleotideColumn_DNA5::INSERT_SEQUENCE
    and variant.mutant() != NucleotideColumn_DNA5::DELETE_NUCLEOTIDE;

  }
  bool applyFilter(const CodonDelete& variant) const override { return true; }

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
  bool applyFilter(const ReadCountVariant& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const SNPVariantDNA5& variant) const override { return not filter_ptr_->applyFilter(variant); }
  bool applyFilter(const CodonDelete& variant) const override { return not filter_ptr_->applyFilter(variant); }

  std::string filterName() const final { return "NOT (" + filter_ptr_->filterName() + ")"; }

private:

  std::unique_ptr<const T> filter_ptr_;

};




}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_FILTER_H
