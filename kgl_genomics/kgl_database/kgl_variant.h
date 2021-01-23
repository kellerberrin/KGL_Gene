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
#include "kgl_genome_collection.h"
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
                  ContigOffset_t contig_offset,
                  bool passedFilters) : contig_id_(std::move(contig_id)),
                                        phase_id_(phase_id),
                                        contig_offset_(contig_offset),
                                        pass_(passedFilters) {}
  virtual ~VariantSequence() = default;

  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  [[nodiscard]] ContigOffset_t offset() const { return contig_offset_; }
  [[nodiscard]] PhaseId_t phaseId() const { return phase_id_; }
  [[nodiscard]] bool passFilters() const { return pass_; }
  void updatePhaseId(PhaseId_t phase_id) { phase_id_ = phase_id; }


  [[nodiscard]] std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

  constexpr static const PhaseId_t UNPHASED = 255;
  constexpr static const PhaseId_t DIPLOID_PHASE_A = 0;  // By convention the female derived contig is first.
  constexpr static const PhaseId_t DIPLOID_PHASE_B = 1;
  constexpr static const PhaseId_t HAPLOID_PHASED = 2;

private:

  const ContigId_t contig_id_;                          // The contig of this variant
  PhaseId_t phase_id_;                                  // The phase of this variant (which homologous contig)
  const ContigOffset_t contig_offset_;                  // Location on the contig.
  bool pass_;                                           // True if the VCF record was annotated with "PASS" filters.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This variant class reflects the variant information presented in VCF files.
// Previous variant objects based on the SAM/BAM file format have been superseded by this class.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { INDEL, TRANSITION, TRANSVERSION };

class Variant : public VariantSequence {

public:

  Variant(   const ContigId_t& contig_id,
             PhaseId_t phase_id,
             ContigOffset_t contig_offset,
             bool passedFilters,
             const VariantEvidence& evidence,
             StringDNA5&& reference,
             StringDNA5&& alternate,
             std::string identifier)
  : VariantSequence(contig_id, phase_id, contig_offset, passedFilters),
    evidence_(evidence),
    reference_(std::move(reference)),
    alternate_(std::move(alternate)),
    identifier_(std::move(identifier)) {}

  ~Variant() override = default;

  // Create a copy of the variant on heap.
  // Important - all the original variant evidence is also attached to the new variant object.
  [[nodiscard]] std::unique_ptr<Variant> clone() const;

  [[nodiscard]] size_t alternateSize() const { return alternate_.length(); }

  [[nodiscard]] size_t referenceSize() const { return reference_.length(); }

  [[nodiscard]] const std::string& identifier() const { return identifier_; }

  [[nodiscard]] VariantType variantType() const;

  [[nodiscard]] bool isSNP() const { return reference().length() == 1 and alternate().length() == 1; }

  [[nodiscard]] bool equivalent(const Variant& cmp_var) const; // Same phase

  [[nodiscard]] bool homozygous(const Variant& cmp_var) const; // Different phase

  [[nodiscard]] bool analogous(const Variant& cmp_var) const; // No phase test

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

  // Set the alternate to the reference; invalidates the evidence structure.
  // Used to specify major alleles (no genome change).
  [[nodiscard]] std::unique_ptr<Variant> cloneNullVariant() const;
  // Checks if the reference and alternate are equivalent.
  [[nodiscard]] bool isNullVariant() const { return reference_ == alternate_; }

private:

  const VariantEvidence evidence_;                      // VCF File based information payload about this variant
  const DNA5SequenceLinear reference_;                  // reference sequence (ref allele)
  const DNA5SequenceLinear alternate_;                  // alternate sequence (alt allele)
  const std::string identifier_;                      // The VCF supplied variant identifier.

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
