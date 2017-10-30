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
#include "kgl_amino.h"
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
class CodonDelete; // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  virtual bool applyFilter(const Variant& variant) const = 0;
  virtual bool applyFilter(const ReadCountVariant& variant) const = 0;
  virtual bool applyFilter(const SNPVariantDNA5& variant) const = 0;
  virtual bool applyFilter(const CodonDelete& variant) const = 0;

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

private:

  virtual bool applyFilter(const VariantFilter& filter) const { return filter.applyFilter(*this); }
  virtual bool equivalent(const Variant& cmp_var) const = 0;

};



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

private:

  const CompoundVariantMap variant_map_;


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A compound delete. 3 contiguous SNP deletions aligned on a Gene codon boundary delete a single Amino Acid
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CodonDelete : public CompoundVariant {

public:

  CodonDelete(std::shared_ptr<const ContigFeatures> contig_ptr,
              ContigOffset_t contig_offset,
              const CompoundVariantMap variant_map,
              ContigOffset_t deleted_codon,
              AminoAcidTypes::AminoType deleted_amino) : CompoundVariant(contig_ptr, contig_offset, variant_map),
                                                         deleted_codon_(deleted_codon),
                                                         deleted_amino_(deleted_amino) {}
  ~CodonDelete() override = default;


private:

  ContigOffset_t deleted_codon_;
  AminoAcidTypes::AminoType deleted_amino_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A read count variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ReadCountVariant : public Variant {

public:

  ReadCountVariant(const std::shared_ptr<const ContigFeatures> contig_ptr,
                   ContigOffset_t contig_offset,
                   NucleotideReadCount_t read_count,
                   NucleotideReadCount_t mutant_count,
                   NucleotideReadCount_t const count_array[],
                   ContigSize_t count_array_size) : Variant(contig_ptr, contig_offset),
                                                            read_count_(read_count),
                                                            mutant_count_(mutant_count) {

    for(ContigOffset_t idx = 0; idx < count_array_size; ++idx) {

      count_array_.push_back(count_array[idx]);

    }

  }
  ReadCountVariant(const ReadCountVariant& variant) = default;
  ~ReadCountVariant() override = default;

  NucleotideReadCount_t readCount() const { return read_count_; }
  NucleotideReadCount_t mutantCount() const { return mutant_count_; }
  const std::vector<NucleotideReadCount_t>& countArray() const { return count_array_; }

  double proportion() const { return static_cast<double>(mutant_count_) / static_cast<double>(read_count_); }

private:

  NucleotideReadCount_t read_count_;
  NucleotideReadCount_t mutant_count_;
  std::vector<NucleotideReadCount_t> count_array_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  A simple SNP variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class SNPVariantDNA5 : public ReadCountVariant {

public:

  SNPVariantDNA5(const std::shared_ptr<const ContigFeatures> contig_ptr,
                 ContigOffset_t contig_offset,
                 NucleotideReadCount_t read_count,
                 NucleotideReadCount_t mutant_count,
                 NucleotideReadCount_t const count_array[],
                 ContigSize_t  count_array_size,
                 typename NucleotideColumn_DNA5::NucleotideType reference,
                 typename NucleotideColumn_DNA5::NucleotideType mutant)
  : ReadCountVariant(contig_ptr, contig_offset, read_count, mutant_count, count_array, count_array_size),
    reference_(reference),
    mutant_(mutant) {}

  SNPVariantDNA5(const SNPVariantDNA5& variant) = default;

  ~SNPVariantDNA5() override = default;

  bool equivalent(const Variant& cmp_var) const override;

  const typename NucleotideColumn_DNA5::NucleotideType& reference() const { return reference_; }
  const typename NucleotideColumn_DNA5::NucleotideType& mutant() const { return mutant_; }

  std::string output() const override;

private:

  typename NucleotideColumn_DNA5::NucleotideType reference_;
  typename NucleotideColumn_DNA5::NucleotideType mutant_;

  bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigVariant - All the variant features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetVariantMap = std::multimap<ContigOffset_t, std::shared_ptr<const Variant>>;
class ContigVariant {

public:

  explicit ContigVariant(const ContigId_t& contig_id) : contig_id_(contig_id) {}
  ContigVariant(const ContigVariant&) = default;
  ~ContigVariant() = default;

  ContigVariant& operator=(const ContigVariant&) = default;

  void addVariant(ContigOffset_t contig_offset, std::shared_ptr<const Variant>& variant_ptr);
  const ContigId_t& contigId() const { return contig_id_; }
  size_t variantCount() const { return offset_variant_map_.size(); }

  // Set functions.
  bool isElement(const Variant& variant) const;
  std::shared_ptr<ContigVariant> Union(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Intersection(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;
  std::shared_ptr<ContigVariant> Difference(std::shared_ptr<const ContigVariant> contig_variant_ptr) const;

  std::shared_ptr<ContigVariant> filterVariants(const VariantFilter& filter) const;

  const OffsetVariantMap& getMap() const { return offset_variant_map_; }

  friend std::ostream& operator<<(std::ostream &os, const ContigVariant& contig_variant);

private:

  ContigId_t contig_id_;
  OffsetVariantMap offset_variant_map_;

};

std::ostream & operator<<(std::ostream &os, const ContigVariant& contig_variant);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeVariant - A map of contig variants
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GenomeVariantMap = std::map<ContigId_t, std::shared_ptr<ContigVariant>>;
class GenomeVariant {

public:

  explicit GenomeVariant(const VariantType_t& variant_type, const GenomeId_t& genome_id) : variant_type_(variant_type),
                                                                                           genome_id_(genome_id) {}
  GenomeVariant(const GenomeVariant&) = default;
  ~GenomeVariant() = default;

  GenomeVariant& operator=(const GenomeVariant& genome_variant) = default;

  const GenomeId_t& genomeId() const { return genome_id_; }
  void genomeId(const GenomeId_t& contig_id) { genome_id_ = contig_id; }

  size_t contigCount() const { return genome_variant_map_.size(); }

  bool addContigVariant(std::shared_ptr<ContigVariant>& contig_variant);
  std::shared_ptr<GenomeVariant> filterVariants(const VariantFilter& filter) const;

  const GenomeVariantMap& contigMap() const { return genome_variant_map_; }

  bool isElement(const Variant& variant) const;
  std::shared_ptr<GenomeVariant> Union(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Intersection(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;
  std::shared_ptr<GenomeVariant> Difference(std::shared_ptr<const GenomeVariant> genome_variant_ptr) const;

  friend std::ostream & operator<<(std::ostream &os, const GenomeVariant& genome_variant);

private:

  VariantType_t variant_type_;
  GenomeId_t genome_id_;
  GenomeVariantMap genome_variant_map_;

};

std::ostream & operator<<(std::ostream &os, const GenomeVariant& genome_variant);


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
