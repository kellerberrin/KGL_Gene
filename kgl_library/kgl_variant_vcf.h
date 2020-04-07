//
// Created by kellerberrin on 31/08/18.
//

#ifndef KGL_VARIANT_VCF_H
#define KGL_VARIANT_VCF_H


#include <map>
#include <memory>
#include <vector>
#include <sstream>
#include "kgl_genome_types.h"
#include "kgl_alphabet_amino.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant.h"
#include "kgl_genome_db.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This variant class reflects the variant information presented in VCF files.
// Previous variant objects based on the SAM/BAM file format have been superseded by this class.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VCFVariant : public Variant {

public:

  VCFVariant(const GenomeId_t& genome_id,
             const ContigId_t& contig_id,
             PhaseId_t phase_id,
             ContigOffset_t contig_offset,
             std::shared_ptr<const VariantEvidence> evidence_ptr,
             StringDNA5&& reference,
             StringDNA5&& alternate)
  : Variant(genome_id, contig_id, phase_id, contig_offset, evidence_ptr),
    reference_(std::move(reference)),
    alternate_(std::move(alternate)) {}

  ~VCFVariant() override = default;

  [[nodiscard]] size_t size() const override { return alternate_.length(); }

  [[nodiscard]] size_t referenceSize() const override { return reference_.length(); }

  [[nodiscard]] VariantType variantType() const override { return VariantType::VCF_VARIANT; }

  [[nodiscard]] bool isSNP() const override { return reference().length() == 1 and alternate().length() == 1; }

  [[nodiscard]] bool equivalent(const Variant& cmp_var) const override;

  [[nodiscard]] bool lessThan(const Variant& cmp_var) const override;

  [[nodiscard]] const DNA5SequenceLinear& reference() const { return reference_; }

  [[nodiscard]] const DNA5SequenceLinear& alternate() const { return alternate_; }

  [[nodiscard]] std::string output(char delimiter, VariantOutputIndex output_index, bool detail) const override;

  [[nodiscard]] std::string mutation(char delimiter, VariantOutputIndex output_index) const override;

  [[nodiscard]] bool mutateSequence( SignedOffset_t offset_adjust,
                                     std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                     SignedOffset_t& sequence_size_modify) const override;


private:

  const DNA5SequenceLinear reference_;
  const DNA5SequenceLinear alternate_;

  [[nodiscard]] bool applyFilter(const VariantFilter& filter) const override { return filter.applyFilter(*this); }

  // Generate a CIGAR by comparing the reference to the alternate.
  [[nodiscard]] std::string alternateCigar() const;

  [[nodiscard]] size_t alternateSize(size_t reference_size) const;

  // Mutate a sequence by adding and subtracting subsequences at the designated offset
  [[nodiscard]] static bool performMutation( ContigOffset_t offset,
                                             std::shared_ptr<DNA5SequenceLinear> mutated_sequence_ptr,
                                             const DNA5SequenceLinear& delete_subsequence,
                                             const DNA5SequenceLinear& add_subsequence);

  [[nodiscard]] bool preceedingMutation( SignedOffset_t adjusted_offset,
                                         std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                         SignedOffset_t& sequence_size_modify) const;

};


}   // end namespace









#endif //KGL_VARIANT_VCF_H
