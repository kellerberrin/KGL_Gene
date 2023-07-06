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
             DNA5SequenceLinear&& reference,
             DNA5SequenceLinear&& alternate,
             const VariantEvidence& evidence) :
      contig_id_(std::move(contig_id)),
      contig_reference_offset_(contig_reference_offset),
      phase_id_(phase_id),
      identifier_(std::move(identifier)),
      reference_(std::move(reference)),
      alternate_(std::move(alternate)),
      evidence_(evidence) { ++object_count_; }

  ~Variant() { --object_count_; }

  // Create a copy of the variant on heap.
  // Important - all the original variant evidence is also attached to the new variant object.
  [[nodiscard]] std::unique_ptr<Variant> clone() const;
  // Set the alternate to the reference; invalidates the evidence structure.
  // Used to specify major alleles (no genome change).
  [[nodiscard]] std::unique_ptr<Variant> cloneNullVariant() const;
  // Clone the variant using canonical reference and alternate.
  // SNPs are represented as '1X', deletes as '1MnD' and inserts as '1MnI'.
  // The first argument is the adjusted reference, the second is the adjusted alternate,
  // the third is the adjusted variant offset. The adjusted variant offset always indicates the offset
  // at which the first nucleotide of the canonical reference and canonical alternate is found.
  [[nodiscard]] std::unique_ptr<Variant> cloneCanonical() const;
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


  [[nodiscard]] bool operator<(const Variant& cmp_var) const { return lessThan(cmp_var); };
  [[nodiscard]] bool operator==(const Variant& cmp_var) const { return equivalent(cmp_var); };

  [[nodiscard]] bool filterVariant(const BaseFilter& filter) const;

  // Location specific parameters.
  // The identifier of contiguous region (chromosome or scaffold) where the variant is located.
  [[nodiscard]] const ContigId_t& contigId() const { return contig_id_; }
  // The offset used for storing the allele in the database.
  // The ZERO BASED offset of the allele in the VCF file. Note, this is NOT the 1 based offset used in VCFs, Gffs etc.
  [[nodiscard]] ContigOffset_t offset() const { return contig_reference_offset_; }
  // For a sequence on the interval [a, b), Given a start offset a and a size (b-a). Determine if the variant will
  // modify the sequence. Note that this is different to just translating the sequence offsets. Any upstream indel will
  // modify the sequence [a, b) offsets but may not actually modify any of the nucleotides in the sequence.
  [[nodiscard]] bool sequenceModifier(ContigOffset_t sequence_start, ContigSize_t sequence_size) const;
  // The extentOffset() of the variant is used to assess if a CANONICAL variant modifies a particular region
  // of a sequence in the interval [a, b). The offset is the canonical offset (see canonicalSequences()) and the extent
  // is 1 for a (canonical) SNP and insert. A delete extent is the number of deleted nucleotides (ref.length() - alt.length() + 1).
  [[nodiscard]] std::pair<ContigOffset_t, ContigSize_t> extentOffset() const;

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
  [[nodiscard]] std::string cigar() const;
  // Reduces the variant to a canonical reference and alternate format.
  // SNPs are represented as '1X', deletes as '1MnD' and inserts as '1MnI'.
  // The first argument is the adjusted reference, the second is the adjusted alternate,
  // the third is the adjusted variant offset. The adjusted variant offset always indicates the offset
  // at which the first nucleotide of the canonical reference and canonical alternate is found.
  [[nodiscard]] std::tuple<DNA5SequenceLinear, DNA5SequenceLinear, ContigOffset_t> canonicalSequences() const;

  // Used to check memory usage and identify any memory leaks.
  [[nodiscard]] static size_t objectCount() { return object_count_; }

private:

  const ContigId_t contig_id_;                          // The contig of this variant
  const ContigOffset_t contig_reference_offset_;        // Physical Location of the start of the reference sequence the contig.
  const VariantPhase phase_id_;                         // The phase of this variant (which homologous contig)
  const std::string identifier_;                        // The VCF supplied variant identifier such as (HSapien) "rs187084".
  const DNA5SequenceLinear reference_;                  // reference sequence (ref allele)
  const DNA5SequenceLinear alternate_;                  // alternate sequence (alt allele)
  const VariantEvidence evidence_;                      // VCF File based information payload about this variant

  inline static std::atomic<size_t> object_count_{0};  // Used to check memory usage and identify any memory leaks.


  // Common (same nucleotide) Prefix and Suffix size for refernce and alternate.
  [[nodiscard]] size_t commonPrefix() const { return reference().commonPrefix(alternate()); }
  [[nodiscard]] size_t commonSuffix() const { return reference().commonSuffix(alternate()); }

};


// Used in the processAll()) templates.
template<class ObjFunc>
using MemberVariantFunc = bool (ObjFunc::*)(const std::shared_ptr<const Variant>&);
using VariantProcessFunc = std::function<bool(const std::shared_ptr<const Variant>&)>;


}   // end namespace



#endif //KGL_VARIANT_H
