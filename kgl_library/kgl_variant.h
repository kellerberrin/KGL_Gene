//
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



namespace kellerberrin::genome {   //  organization level namespace


class Variant;   // Forward decl.

class VariantFilter {

public:

  VariantFilter() = default;
  virtual ~VariantFilter() = default;

  [[nodiscard]] virtual bool applyFilter(const Variant& variant) const = 0;

  [[nodiscard]] virtual std::shared_ptr<VariantFilter> clone() const = 0;

  [[nodiscard]] std::string filterName() const { return filter_name_; }

  void filterName(std::string filter_name) { filter_name_ = std::move(filter_name); }


private:

  std::string filter_name_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class VariantSequenceType { CDS_CODING, INTRON, NON_CODING };
class VariantSequence {

public:

  VariantSequence(ContigId_t contig_id,
                  PhaseId_t phase_id,
                  ContigOffset_t contig_offset) : contig_id_(std::move(contig_id)),
                                                  phase_id_(phase_id),
                                                  contig_offset_(contig_offset) {}
  virtual ~VariantSequence() = default;

  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] ContigOffset_t offset() const { return contig_offset_; }
  [[nodiscard]] PhaseId_t phaseId() const { return phase_id_; }
  virtual void updatePhaseId(PhaseId_t phase_id) { protectedPhaseId(phase_id); }


  [[nodiscard]] std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

  constexpr static const PhaseId_t UNPHASED = 255;
  constexpr static const PhaseId_t DIPLOID_PHASE_F = 0;  // By convention the female derived contig is first.
  constexpr static const PhaseId_t DIPLOID_PHASE_M = 1;

protected:

  void protectedPhaseId(PhaseId_t phase_id) { phase_id_ = phase_id; }

private:

  const ContigId_t contig_id_;                          // The contig of this variant
  PhaseId_t phase_id_;                                  // The phase of this variant (which homologous contig)
  const ContigOffset_t contig_offset_;                  // Location on the contig.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This variant class reflects the variant information presented in VCF files.
// Previous variant objects based on the SAM/BAM file format have been superseded by this class.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { VCF_VARIANT };

class Variant : public VariantSequence {

public:

  Variant(   const ContigId_t& contig_id,
             PhaseId_t phase_id,
             ContigOffset_t contig_offset,
             const VariantEvidence& evidence,
             StringDNA5&& reference,
             StringDNA5&& alternate)
  : VariantSequence(contig_id, phase_id, contig_offset),
    evidence_(evidence),
    reference_(std::move(reference)),
    alternate_(std::move(alternate)) {}

  ~Variant() override = default;

  [[nodiscard]] size_t alternateSize() const { return alternate_.length(); }

  [[nodiscard]] size_t referenceSize() const { return reference_.length(); }

  [[nodiscard]] VariantType variantType() const { return VariantType::VCF_VARIANT; }

  [[nodiscard]] bool isSNP() const { return reference().length() == 1 and alternate().length() == 1; }

  [[nodiscard]] bool equivalent(const Variant& cmp_var) const;

  [[nodiscard]] bool lessThan(const Variant& cmp_var) const;

  [[nodiscard]] const DNA5SequenceLinear& reference() const { return reference_; }

  [[nodiscard]] const DNA5SequenceLinear& alternate() const { return alternate_; }

  [[nodiscard]] const VariantEvidence& evidence() const { return evidence_; }

  [[nodiscard]] std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const;

  [[nodiscard]] std::string mutation(char delimiter, VariantOutputIndex output_index) const;

  [[nodiscard]] bool mutateSequence( SignedOffset_t offset_adjust,
                                     DNA5SequenceLinear& dna_sequence,
                                     SignedOffset_t& sequence_size_modify) const;

  [[nodiscard]] bool operator<(const Variant& cmp_var) const { return lessThan(cmp_var); };
  [[nodiscard]] bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

  [[nodiscard]] std::string name() const;

  [[nodiscard]] bool filterVariant(const VariantFilter& filter) const { return filter.applyFilter(*this); }

private:

  const VariantEvidence evidence_;                      // VCF File based information payload about this variant
  const DNA5SequenceLinear reference_;                  // reference sequence (allele)
  const DNA5SequenceLinear alternate_;                  // alternate sequence (allele)

  // Generate a CIGAR by comparing the reference to the alternate.
  [[nodiscard]] std::string alternateCigar() const;

  [[nodiscard]] size_t alternateSize(size_t reference_size) const;

  // Mutate a sequence by adding and subtracting subsequences at the designated offset
  [[nodiscard]] static bool performMutation( ContigOffset_t offset,
                                             DNA5SequenceLinear& mutated_sequence,
                                             const DNA5SequenceLinear& delete_subsequence,
                                             const DNA5SequenceLinear& add_subsequence);

  [[nodiscard]] bool preceedingMutation( SignedOffset_t adjusted_offset,
                                         DNA5SequenceLinear& dna_sequence,
                                         SignedOffset_t& sequence_size_modify) const;

};


}   // end namespace



#endif //KGL_VARIANT_H
