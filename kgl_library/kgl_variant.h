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
#include "kgl_variant_evidence.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The abstract VariantFilter class uses the visitor pattern.
// Concrete variant filters are defined in kgl_filter.h
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class VCFVariant;   // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  virtual bool applyFilter(const VCFVariant& variant) const = 0;

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
                  std::shared_ptr<const VariantEvidence> evidence_ptr) : genome_id_(genome_id),
                                                                         contig_id_(contig_id),
                                                                         phase_id_(phase_id),
                                                                         contig_offset_(contig_offset),
                                                                         evidence_ptr_(evidence_ptr) {}
  virtual ~VariantSequence() = default;

  const GenomeId_t& genomeId() const { return genome_id_; }
  const ContigId_t& contigId() const { return contig_id_; }
  ContigOffset_t offset() const { return contig_offset_; }
  PhaseId_t phaseId() const { return phase_id_; }
  virtual void updatePhaseId(PhaseId_t phase_id) { protectedPhaseId(phase_id); }

  std::shared_ptr<const VariantEvidence> evidence() const { return evidence_ptr_; }

  std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

  static constexpr PhaseId_t UNPHASED = 255;

protected:

  void protectedPhaseId(PhaseId_t phase_id) { phase_id_ = phase_id; }

private:

  GenomeId_t genome_id_;                          // The source of this variant
  ContigId_t contig_id_;                          // The contig of this variant
  PhaseId_t phase_id_;                            // The phase of this variant (which homologous contig)
  ContigOffset_t contig_offset_;                  // Location on the contig.
  const std::shared_ptr<const VariantEvidence> evidence_ptr_;  // Addition information about the variant.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Base class for a genome variant. Modelled on the VCF file format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { VCF_VARIANT };


class Variant : public VariantSequence {

public:


  Variant(const GenomeId_t& genome_id,
          const ContigId_t& contig_id,
          PhaseId_t phase_id,
          ContigOffset_t contig_offset,
          std::shared_ptr<const VariantEvidence> evidence_ptr) : VariantSequence(genome_id,
                                                                                 contig_id,
                                                                                 phase_id,
                                                                                 contig_offset,
                                                                                 evidence_ptr) {}
  ~Variant() override = default;
  bool filterVariant(const VariantFilter& filter) const { return applyFilter(filter); }

  std::string name() const;

  virtual VariantType variantType() const = 0;
  virtual size_t size() const = 0;
  virtual size_t referenceSize() const { return 1; }
  virtual std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const = 0;
  virtual std::string mutation(char delimiter, VariantOutputIndex output_index) const = 0;
  virtual bool mutateSequence(SignedOffset_t offset_adjust,
                              std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                              SignedOffset_t& sequence_size_modify) const = 0;

  bool isCompound() const { return size() > 1; }
  bool isSingle() const { return size() == 1; }
  virtual bool isSNP() const { return false; }

  virtual bool equivalent(const Variant& cmp_var) const = 0;
  bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

  virtual bool lessThan(const Variant& cmp_var) const = 0;
  bool operator<(const Variant& cmp_var) const { return lessThan(cmp_var); };

  virtual std::shared_ptr<Variant> clone() const = 0;  // Polymorphic copy constructor

private:

  virtual bool applyFilter(const VariantFilter& filter) const = 0;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_H
