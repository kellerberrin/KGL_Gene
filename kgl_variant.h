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
#include "kgl_sequence_base.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class SNPVariant; // Forward decl.
class DeleteVariant; // Forward decl.
class InsertVariant; // Forward decl.
class CompoundDelete; // Forward decl.
class CompoundInsert; // Forward decl.
class CompoundSNP; // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  virtual bool applyFilter(const SNPVariant& variant) const = 0;
  virtual bool applyFilter(const DeleteVariant& variant) const = 0;
  virtual bool applyFilter(const InsertVariant& variant) const = 0;
  virtual bool applyFilter(const CompoundDelete& variant) const = 0;
  virtual bool applyFilter(const CompoundInsert& variant) const = 0;
  virtual bool applyFilter(const CompoundSNP& variant) const = 0;

  virtual std::string filterName() const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant offset output convention.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used to display output as 1 based or zero based offsets. Never, ever, use START_1_BASED internally.
enum class VariantOutputIndex { START_1_BASED, START_0_BASED};   // Used for output functions - default START_1_BASED

// helper function - only ever use for output.
std::string offsetOutput(ContigOffset_t offset, VariantOutputIndex output_base);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class VariantSequenceType { CDS_CODING, INTRON, NON_CODING };
class VariantSequence {

public:

  VariantSequence(const std::string& variant_source,
                  std::shared_ptr<const ContigFeatures>  contig_ptr,
                  ContigOffset_t contig_offset,
                  Phred_t quality) : variant_source_(variant_source),
                                     contig_ptr_(contig_ptr),
                                     contig_offset_(contig_offset),
                                     quality_(quality) { }
  virtual ~VariantSequence() = default;

  const std::string& variantSource() const { return variant_source_; }

  std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

  std::shared_ptr<const ContigFeatures> contig() const { return contig_ptr_; }
  ContigOffset_t offset() const { return contig_offset_; }

  VariantSequenceType type() const;
  std::string typestr() const;

  Phred_t quality() const { return quality_; }
  void quality(Phred_t quality_update) { quality_ = quality_update; }

  FeatureIdent_t codingSequenceId() const {

    return type() == VariantSequenceType::CDS_CODING ? codingSequences().getFirst()->getCDSParent()->id() : NULL_ID;

  }

  const CodingSequenceArray& codingSequences() const { return coding_sequences_; }
  const GeneVector& geneMembership() const { return gene_membership_; }

  // returns false if not in a coding sequence.
  bool codonOffset(ContigOffset_t& codon_offset, ContigSize_t& base_in_codon) const;

  void defineIntron(std::shared_ptr<const GeneFeature> gene_ptr);
  void defineCoding(std::shared_ptr<const CodingSequence> coding_sequence_ptr);
  void defineNonCoding();

  static constexpr const char* NULL_ID = "NULL";

private:

  std::string variant_source_;                          // The source of this variant
  std::shared_ptr<const ContigFeatures> contig_ptr_;    // The contig.
  ContigOffset_t contig_offset_;                        // Location on the contig.
  GeneVector gene_membership_;                          // Membership includes introns (empty for non-coding)
  CodingSequenceArray coding_sequences_;  // Coding sequence for variant (empty for introns and non-coding)
  Phred_t quality_;                       // Phred (-10log10) quality that the variant is not valid.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { SNP, DELETE, INSERT, COMPOUND_SNP, COMPOUND_INSERT, COMPOUND_DELETE };

class Variant : public VariantSequence {

public:


  Variant(const std::string& variant_source,
          const std::shared_ptr<const ContigFeatures> contig_ptr,
          ContigOffset_t contig_offset,
          Phred_t quality) : VariantSequence(variant_source, contig_ptr, contig_offset, quality) {}
  ~Variant() override = default;
  bool filterVariant(const VariantFilter& filter) const { return applyFilter(filter); }


  const ContigId_t& contigId() const { return contig()->contigId(); }
  ContigOffset_t contigOffset() const { return offset(); }
  bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };
  bool offsetOverlap(const Variant& cmp_var) const;  // Particuarly relevant for compound variants.

  std::string name() const;

  virtual VariantType variantType() const = 0;
  virtual size_t size() const = 0;
  virtual std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const = 0;
  virtual std::string mutation(char delimiter, VariantOutputIndex output_index) const = 0;
  virtual bool mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                    SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                                    ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                                    SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                                    std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const = 0;

  bool isCompound() const { return size() > 1; }
  bool isSingle() const { return size() == 1; }
  bool isSNP() const { return variantType() == VariantType::SNP or variantType() == VariantType::COMPOUND_SNP; }
  bool isDelete() const { return variantType() == VariantType::DELETE or variantType() == VariantType::COMPOUND_DELETE; }
  bool isInsert() const { return variantType() == VariantType::INSERT or variantType() == VariantType::COMPOUND_INSERT; }

  virtual bool equivalent(const Variant& cmp_var) const = 0;
  virtual std::shared_ptr<Variant> clone() const = 0;  // Polymorphic copy constructor

protected:

  static constexpr const char CODON_BASE_SEPARATOR = ':';

private:

  virtual bool applyFilter(const VariantFilter& filter) const = 0;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
