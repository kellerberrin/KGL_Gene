//
// Created by kellerberrin on 13/10/17.
//

#ifndef KGL_VARIANT_H
#define KGL_VARIANT_H

#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include <algorithm>
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

  [[nodiscard]] bool mutateSequence( ContigOffset_t allele_offset,
                                     SignedOffset_t offset_adjust,
                                     DNA5SequenceLinear& dna_sequence,
                                     SignedOffset_t& sequence_size_modify) const;

  [[nodiscard]] bool operator<(const Variant& cmp_var) const { return lessThan(cmp_var); };
  [[nodiscard]] bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

  [[nodiscard]] bool filterVariant(const BaseFilter& filter) const;

  // Location specific parameters.
  // The identifier of contiguous region (chromosome or scaffold) where the variant is located.
  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  // The offset used for storing the allele in the database.
  // The ZERO BASED offset of the allele in the VCF file. Note, this is NOT the 1 based offset used in VCFs, Gffs etc.
  [[nodiscard]] ContigOffset_t offset() const { return contig_reference_offset_; }
  // Actual offset of the allele. For SNPs this will (generally but not always) be zero. For indels this will (generally) be >= 1.
  // The offset is the number of identical nucleotides at the beginning of both the reference and alternate.
  // For example; if we have a variant cigars of ("5M4D", "5M6I" or "5M1X"), the allele offset will be 5 in all cases.
  [[nodiscard]] AlleleOffset_t alleleOffset() const { return contig_allele_offset_; }
  // For a sequence on the interval [a, b), Given a start offset a and a size (b-a). Determine if the variant will
  // modify the sequence. Note that this is different to just translating the sequence offsets. Any upstream indel will
  // modify the sequence [a, b) offsets but may not actually modify any of the nucleotides in the sequence.
  [[nodiscard]] bool sequenceModifier(ContigOffset_t sequence_start, ContigSize_t sequence_size) const;
  // The extentOffset() of the variant is used to assess if a variant modifies a particular region of a sequence in the interval [a, b).
  // With SNP variants the extentOffset().first is c = (offset() + alleleOffset()) so for "5M1X" the extent offset will be c = (offset() + 5).
  // The extent size (extentOffset().second) will be 1. Thus extentOffset() will return the pair [offset()+5,1].
  // The delete variant "5M4D" will have an extentOffset().first of (offset() + 5) and an extent size of 4 (reference.length() - 5).
  // Thus the delete variant "5M4D" can modify the sequence [a, b) if it is outside (before) the interval, if for example, (offset() + alleleOffset()) = a-2.
  // The deleted nucleotides will be {a-2. a-1, a,  a+1}. The delete variant can preceed the interval [a, b).
  // The insert variant "5M6I", the same as an SNP, will have an extentOffset() of 5 and an extent size of 1.
  // Unlike a delete variant, the insert variant will only modify [a, b) if (offset() + alleleOffset()) is in [a, b).
  // All variants will modify a sequence [a, b) (not just translate it's offsets) if the following condition is met:
  // bool modified = (extentOffset().first + extentOffset.second) > a or (extentOffset().first < b);
  [[nodiscard]] std::pair<ContigOffset_t, ContigSize_t> extentOffset() const;
  // The actual start of the allele modification of the underlying contig sequence, is >= 0.
  [[nodiscard]] ContigOffset_t alleleMutateOffset() const { return offset()  + alleleOffset(); }

  [[nodiscard]] VariantPhase phaseId() const { return phase_id_; }
  // An assigned variant reference such as (HSapien) "rs187084".
  [[nodiscard]] const std::string& identifier() const { return identifier_; }

  // A string hash unique upto phase (not phase specific). In HGVS format.
  // This string hash is extensively used as a map key for variants.
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
  // Prefix and Suffix size
  [[nodiscard]] size_t commonPrefix() const { return reference().commonPrefix(alternate()); }
  [[nodiscard]] size_t commonSuffix() const { return reference().commonSuffix(alternate()); }
  // Remove common prefix and suffix nucleotides from the reference and alternate.
  // .first is the trimmed reference, .second is the trimmed alternate.
  [[nodiscard]] std::pair<DNA5SequenceLinear, DNA5SequenceLinear> trimmedSequences() const;

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

  inline static std::atomic<size_t> object_count_{0};  // Used to check memory usage and identify any memory leaks.


  [[nodiscard]] size_t alternateSize(size_t reference_size) const;


  // Mutate a sequence by adding and subtracting subsequences at the designated offset
  [[nodiscard]] static bool performMutation( ContigOffset_t offset,
                                             DNA5SequenceLinear& mutated_sequence,
                                             const DNA5SequenceLinear& delete_subsequence,
                                             const DNA5SequenceLinear& add_subsequence,
                                             SignedOffset_t& sequence_size_modify);

  [[nodiscard]] bool preceedingMutation( SignedOffset_t adjusted_offset,
                                         DNA5SequenceLinear& dna_sequence,
                                         SignedOffset_t& sequence_size_modify) const;



};


// Used in the processAll()) templates.
template<class ObjFunc>
using MemberVariantFunc = bool (ObjFunc::*)(const std::shared_ptr<const Variant>&);
using VariantProcessFunc = std::function<bool(const std::shared_ptr<const Variant>&)>;


}   // end namespace



#endif //KGL_VARIANT_H
