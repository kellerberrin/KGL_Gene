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

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  virtual bool applyFilter(const SNPVariant& variant) const = 0;
  virtual bool applyFilter(const DeleteVariant& variant) const = 0;
  virtual bool applyFilter(const InsertVariant& variant) const = 0;
  virtual bool applyFilter(const CompoundDelete& variant) const = 0;
  virtual bool applyFilter(const CompoundInsert& variant) const = 0;

  virtual std::string filterName() const = 0;
  virtual std::shared_ptr<VariantFilter> clone() const = 0;

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant offset output convention.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class VariantSequenceType { CDS_CODING, INTRON, NON_CODING };
class VariantSequence {

public:

  VariantSequence(const GenomeId_t& genome_id,
                  const ContigId_t& contig_id,
                  PhaseId_t phase_id,
                  ContigOffset_t contig_offset,
                  Phred_t quality) : genome_id_(genome_id),
                                     contig_id_(contig_id),
                                     phase_id_(phase_id),
                                     contig_offset_(contig_offset),
                                     quality_(quality) { }
  virtual ~VariantSequence() = default;

  const GenomeId_t& genomeId() const { return genome_id_; }
  const ContigId_t& contigId() const { return contig_id_; }
  PhaseId_t phaseId() const { return phase_id_; }
  void phaseId(PhaseId_t phase_id) { phase_id_ = phase_id; }

  std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

  ContigOffset_t offset() const { return contig_offset_; }

  Phred_t quality() const { return quality_; }
  void quality(Phred_t quality_update) { quality_ = quality_update; }

  static constexpr PhaseId_t UNPHASED = 255;

private:

  GenomeId_t genome_id_;                          // The source of this variant
  ContigId_t contig_id_;                          // The contig of this variant
  PhaseId_t phase_id_;                            // The phase of this variant (which homologous contig)
  ContigOffset_t contig_offset_;                  // Location on the contig.
  Phred_t quality_;                       // Phred (-10log10) quality that the variant is not valid.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { SNP, DELETE, INSERT, COMPOUND_INSERT, COMPOUND_DELETE };


class Variant : public VariantSequence {

public:


  Variant(const GenomeId_t& genome_id,
          const ContigId_t& contig_id,
          PhaseId_t phase_id,
          ContigOffset_t contig_offset,
          Phred_t quality) : VariantSequence(genome_id, contig_id, phase_id, contig_offset, quality) {}
  ~Variant() override = default;
  bool filterVariant(const VariantFilter& filter) const { return applyFilter(filter); }

  bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

  bool offsetOverlap(const Variant& cmp_var) const;  // Particuarly relevant for compound variants.

  std::string name() const;

  virtual VariantType variantType() const = 0;
  virtual size_t size() const = 0;
  virtual std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const = 0;
  virtual std::string mutation(char delimiter, VariantOutputIndex output_index) const = 0;
  virtual bool mutateSequence(SignedOffset_t offset_adjust,
                              std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                              SignedOffset_t& sequence_size_modify) const = 0;

  bool isCompound() const { return size() > 1; }
  bool isSingle() const { return size() == 1; }
  bool isSNP() const { return variantType() == VariantType::SNP;  }
  bool isDelete() const { return variantType() == VariantType::DELETE or variantType() == VariantType::COMPOUND_DELETE; }
  bool isInsert() const { return variantType() == VariantType::INSERT or variantType() == VariantType::COMPOUND_INSERT; }

  virtual bool equivalent(const Variant& cmp_var) const = 0;
  virtual std::shared_ptr<Variant> clone() const = 0;  // Polymorphic copy constructor

private:

  virtual bool applyFilter(const VariantFilter& filter) const = 0;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
