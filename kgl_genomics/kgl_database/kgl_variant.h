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
#include "kgl_runtime_resource.h"
#include "kgl_sequence_base.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_filter_virtual.h"
#include "kgl_variant_factory_vcf_parse_cigar.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Variant Chromosome phase
enum class VariantPhase : std::uint8_t { HAPLOID_PHASED = 0,
                                         DIPLOID_PHASE_A = 1, // By convention, the female derived contig is first.
                                         DIPLOID_PHASE_B = 2,
                                         UNPHASED = 255 };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This variant class reflects the variant information presented in VCF files.
// The object is completely constant (cannot be modified). This is because the object is simultaneously
// held by multiple genomes to conserve memory. Therefore modifying a variant would modify the variant everywhere.
// If a modified copy of the variant is required then it must be cloned with modified arguments.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Defined Variant Types.
enum class VariantType { INDEL_DELETE, INDEL_INSERT, TRANSITION, TRANSVERSION };

// Type of variant equality, phased also checks the variant phase.
enum class VariantEquality { PHASED, UNPHASED };

class Variant {

public:

  Variant(   ContigId_t contig_id,
             ContigOffset_t contig_reference_offset,
             VariantPhase phase_id,
             std::string identifier,
             StringDNA5&& reference,
             StringDNA5&& alternate,
             const VariantEvidence& evidence) :
      reference_(std::move(reference)),
      alternate_(std::move(alternate)),
      evidence_(evidence),
      contig_id_(std::move(contig_id)),
      contig_reference_offset_(contig_reference_offset),
      contig_allele_offset_(reference_.commonPrefix(alternate_)),
      phase_id_(phase_id),
      identifier_(std::move(identifier)) { ++object_count_; }

  ~Variant() { --object_count_; }

  // Create a copy of the variant on heap.
  // Important - all the original variant evidence is also attached to the new variant object.
  [[nodiscard]] std::unique_ptr<Variant> clone() const;

  // Set the alternate to the reference; invalidates the evidence structure.
  // Used to specify major alleles (no genome change).
  [[nodiscard]] std::unique_ptr<Variant> cloneNullVariant() const;

  // Clone the variant with a different phase.
  [[nodiscard]] std::unique_ptr<Variant> clonePhase(VariantPhase phase_id) const;

  [[nodiscard]] size_t alternateSize() const { return alternate_.length(); }

  [[nodiscard]] size_t referenceSize() const { return reference_.length(); }

  [[nodiscard]] VariantType variantType() const;

  // Includes ref, alt comparison with 1'X' and n'M' cigars.
  [[nodiscard]] bool isSNP() const;

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

  [[nodiscard]] std::string typeText() const;

  [[nodiscard]] bool filterVariant(const BaseFilter& filter) const;

  // Location specific parameters.
  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  // The offset used for storing the allele in the database. Currently the same as the reference offset.
  [[nodiscard]] ContigOffset_t offset() const { return contig_reference_offset_; }
  // The zero based offset of the allele in the VCF file. Note, this is NOT the 1 based offset used in VCFs, Gffs etc.
  [[nodiscard]] ContigOffset_t referenceOffset() const { return contig_reference_offset_; }
  // Actual offset of the allele. For SNPs this will (generally) be zero. For indels this will (generally) be >= 1.
  [[nodiscard]] AlleleOffset_t alleleOffset() const { return contig_allele_offset_; }
  [[nodiscard]] VariantPhase phaseId() const { return phase_id_; }
  [[nodiscard]] const std::string& identifier() const { return identifier_; }

  // A string hash unique upto phase (not phase specific). In HGVS format.
  [[nodiscard]] std::string HGVS() const;
  // Phase specific hash. In HGVS format plus the variant phase (if applicable - may be unphased) in the format ":B".
  [[nodiscard]] std::string HGVS_Phase() const;

  // The following are comparison functions for ordering variants.
  // Equality hash.
  [[nodiscard]] std::string equalityHash(VariantEquality type) const { return (type == VariantEquality::PHASED) ? HGVS_Phase() : HGVS(); }

  [[nodiscard]] bool equality(const Variant& cmp, VariantEquality type) const { return equalityHash(type) == cmp.equalityHash(type); }

  [[nodiscard]] bool equivalent(const Variant& cmp_var) const { return equality(cmp_var, VariantEquality::PHASED); }

  [[nodiscard]] bool homozygous(const Variant& cmp_var) const { return analogous(cmp_var) and phaseId() != cmp_var.phaseId(); }

  [[nodiscard]] bool analogous(const Variant& cmp_var) const { return equality(cmp_var, VariantEquality::UNPHASED); }

  [[nodiscard]] bool lessThan(const Variant& cmp_var) const {   return HGVS_Phase() < cmp_var.HGVS_Phase(); }

  // Generate a CIGAR by comparing the reference to the alternate.
  [[nodiscard]] std::string alternateCigar() const;

  // Used to check memory usage and identify any memory leaks.
  [[nodiscard]] static size_t objectCount() { return object_count_; }

private:

  const DNA5SequenceLinear reference_;                  // reference sequence (ref allele)
  const DNA5SequenceLinear alternate_;                  // alternate sequence (alt allele)
  const VariantEvidence evidence_;                      // VCF File based information payload about this variant
  const ContigId_t contig_id_;                          // The contig of this variant
  const ContigOffset_t contig_reference_offset_;        // Physical Location of the start of the reference sequence the contig.
  const AlleleOffset_t contig_allele_offset_;           // Offset to where the allele actually occurs. Always 0 for an SNP, always > 0 for an indel.
  const VariantPhase phase_id_;                         // The phase of this variant (which homologous contig)
  const std::string identifier_;                        // The VCF supplied variant identifier.

  inline static std::atomic<size_t> object_count_{0};


  [[nodiscard]] size_t alternateSize(size_t reference_size) const;

  // Mutate a sequence by adding and subtracting subsequences at the designated offset
  [[nodiscard]] static bool performMutation( ContigOffset_t offset,
                                             DNA5SequenceLinear& mutated_sequence,
                                             const DNA5SequenceLinear& delete_subsequence,
                                             const DNA5SequenceLinear& add_subsequence);

  [[nodiscard]] bool preceedingMutation( SignedOffset_t adjusted_offset,
                                         DNA5SequenceLinear& dna_sequence,
                                         SignedOffset_t& sequence_size_modify) const;

  [[nodiscard]] std::string genomeOutput(char delimiter, VariantOutputIndex output_index) const;  // Genome information text.

};


// Used in the processAll()) templates.
template<class ObjFunc>
using MemberVariantFunc = bool (ObjFunc::*)(const std::shared_ptr<const Variant>&);
using VariantProcessFunc = std::function<bool(const std::shared_ptr<const Variant>&)>;


}   // end namespace



#endif //KGL_VARIANT_H
