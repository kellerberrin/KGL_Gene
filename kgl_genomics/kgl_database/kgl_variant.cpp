//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kgl_variant_filter_type.h"
#include "kel_patterns.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::Variant> kgl::Variant::clone() const {

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId(),
                                                                  identifier(),
                                                                  DNA5SequenceLinear(StringDNA5(reference().getAlphabetString())),
                                                                  DNA5SequenceLinear(StringDNA5(alternate().getAlphabetString())),
                                                                  evidence()));

  return variant_ptr;

}


std::unique_ptr<kgl::Variant> kgl::Variant::cloneNullVariant() const {

  VariantEvidence null_evidence; // no evidence is passed through.

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId(),
                                                                  identifier(),
                                                                  DNA5SequenceLinear(StringDNA5(reference().getAlphabetString())),
                                                                  DNA5SequenceLinear(StringDNA5(reference().getAlphabetString())),
                                                                  null_evidence));

  return variant_ptr;

}

std::unique_ptr<kgl::Variant> kgl::Variant::cloneCanonical() const {

  auto [canonical_ref, canonical_alt, canonical_offset] = canonicalSequences();

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  canonical_offset,
                                                                  phaseId(),
                                                                  identifier(),
                                                                  std::move(canonical_ref),
                                                                  std::move(canonical_alt),
                                                                  evidence()));

  return variant_ptr;

}

// Clone with modified phase.
std::unique_ptr<kgl::Variant> kgl::Variant::clonePhase(VariantPhase phaseId) const {

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId,
                                                                  identifier(),
                                                                  DNA5SequenceLinear(StringDNA5(reference().getAlphabetString())),
                                                                  DNA5SequenceLinear(StringDNA5(alternate().getAlphabetString())),
                                                                  evidence()));

  return variant_ptr;

}


bool kgl::Variant::filterVariant(const BaseFilter& filter) const {

  std::shared_ptr<const FilterVariants> variant_filter = std::dynamic_pointer_cast<const FilterVariants>(filter.clone());
  if (variant_filter) {

    return variant_filter->applyFilter(*this);

  }

  ExecEnv::log().error("Variant::filterVariant; Filter: {} is not a variant filter", filter.filterName());
  return false;

}


kgl::VariantType kgl::Variant::variantType() const {

  if (not isSNP()) {

    return (reference().length() <= alternate().length()) ?  VariantType::INDEL_INSERT : VariantType::INDEL_DELETE;

  } else {

    if (DNA5::isTransition(reference().at(0), alternate().at(0))) {

      return VariantType::TRANSITION;

    } else {

     return VariantType::TRANSVERSION;

    }

  }

}


// VCF files can specify SNPs as a cigar such as '4M1X8M'.
// We extend the logic for this possibility.
bool kgl::Variant::isSNP() const {


    // Obvious SNP.
  if (reference_.length() == 1 and alternate_.length() == 1) {

    return true;

  }

  // Not SNP if different lengths.
  if (reference_.length() != alternate_.length()) {

    return false;

  }

  // Check only 1 nucleotide difference.
  bool diff_found{false};
  for (size_t i = 0; i < reference_.length(); ++i) {

    if (reference_[i] != alternate_[i]) {

      if (diff_found) {

        return false;

      }

      diff_found = true;

    }

  }

  return true;

}

// Reduces the variant to a canonical reference and alternate format.
// SNPs are represented as '1X', deletes as '1MnD' and inserts as '1MnI'.
// The first argument is the adjusted reference, the second is the adjusted alternate,
// the third is the adjusted variant offset. The adjusted variant offset always indicates the offset
// at which the first nucleotide of the canonical reference and canonical alternate is found.
std::tuple<kgl::DNA5SequenceLinear, kgl::DNA5SequenceLinear, kgl::ContigOffset_t> kgl::Variant::canonicalSequences() const {

  size_t prefix_size = commonPrefix();
  prefix_size = prefix_size > 0 ? (prefix_size - 1) : 0;  // Adjust for '1MnD' and '1MnI'.
  size_t suffix_size = commonSuffix();
  size_t min_size = std::min(referenceSize(), alternateSize());
  int64_t adj_suffix_size = std::min(min_size - prefix_size - 1, suffix_size); // Adjust for '1MnD' and '1MnI'.
  adj_suffix_size = adj_suffix_size < 0 ? 0 : adj_suffix_size;

  auto canonical_reference = reference().removePrefixSuffix(prefix_size, adj_suffix_size);
  auto canonical_alternate = alternate().removePrefixSuffix(prefix_size, adj_suffix_size);
  ContigOffset_t canonical_offset = offset() + prefix_size;

  return { DNA5SequenceLinear(std::move(canonical_reference)), DNA5SequenceLinear(std::move(canonical_alternate)), canonical_offset };

}

// Check if the variants reference() and alternate() are in canonical form.
// The variant is canonical if SNPs are represented as '1X', deletes as '1MnD' and inserts as '1MnI'.
// A canonical SNP will have referenceSize() == alternateSize() == 1.
// A canonical delete will have referenceSize() > alternateSize() == 1.
// A canonical insert will have 1 == referenceSize() < alternateSize().
bool kgl::Variant::isCanonical() const {

  if (referenceSize() == 1 and referenceSize() == alternateSize()) {
    // Canonical SNP
    return true;

  } else if (alternateSize() == 1 and referenceSize() > alternateSize()) {
    // Canonical Delete
    return true;

  } else if (referenceSize() == 1 and referenceSize() < alternateSize()) {
    // Canonical Insert
    return true;

  }

  return false;

}




// The extentOffset() of the variant is used to assess if a CANONICAL variant modifies a particular region
// of a sequence in the interval [a, b). The offset is the canonical offset (see canonicalSequences()) and the extent
// is 1 for a (canonical) SNP and insert. A delete extent is the number of deleted nucleotides, ref.length().
std::pair<kgl::ContigOffset_t, kgl::ContigSize_t> kgl::Variant::extentOffset() const {

  auto const [canonical_ref, canonical_alt, canonical_offset] = canonicalSequences();
  if (isSNP()) {

    return { canonical_offset, 1 };

  } else {

    if (canonical_ref.length() < canonical_alt.length()) {
      // An insert
      return { canonical_offset, 1 };

    } else if (canonical_ref.length() > canonical_alt.length()) {
      // A delete
      return { canonical_offset, canonical_ref.length()};

    } else {

      ExecEnv::log().warn("Variant::extentOffset; Unexpected variant: {}, cigar:{}, canonical ref: {}, canonical alt: {}, canonical offset: {}"
                          , HGVS(), cigar(), canonical_ref.getSequenceAsString(), canonical_alt.getSequenceAsString(), canonical_offset);

    }

  }

  return {canonical_offset, 1 };

}


// For a sequence on the interval [a, b), Given a start offset a and a size (b-a). Determine if the variant will
// modify the sequence. Note that this is different to just translating the sequence offsets. Any upstream indel will
// modify the sequence [a, b) offsets but may not actually modify any of the nucleotides in the sequence.
bool kgl::Variant::sequenceModifier(ContigOffset_t sequence_start, ContigSize_t sequence_size) const {

  auto [extent_offset, extent_size] = extentOffset();
  return sequence_start < (extent_offset + extent_size) and extent_offset < (sequence_start + sequence_size);

}

std::string kgl::Variant::cigar() const {

  return ParseVCFCigar::generateCigar(reference().getSequenceAsString(), alternate().getSequenceAsString());

}



// Unique upto phase.
std::string kgl::Variant::HGVS() const {

  // 1 is added to the offset to make it 1-based as is the standard in Gffs etc.
  return contigId() + ":g."  + std::to_string(offset() + 1) + reference().getSequenceAsString() + ">" + alternate().getSequenceAsString();

}

// Phase specific hash
std::string kgl::Variant::HGVS_Phase() const {

  // 1 is added to the offset to make it 1-based as is the standard in Gffs etc.
  return contigId() + ":g."  + std::to_string(offset() + 1) + reference().getSequenceAsString() +
         ">" + alternate().getSequenceAsString() + ":" + std::to_string(static_cast<uint8_t>(phaseId()));

}
