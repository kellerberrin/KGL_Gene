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

  StringDNA5 reference_str(reference().getAlphabetString());
  StringDNA5 alternate_str(alternate().getAlphabetString());
  std::string variant_ident = identifier();

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId(),
                                                                  std::move(variant_ident),
                                                                  std::move(reference_str),
                                                                  std::move(alternate_str),
                                                                  evidence()));

  return variant_ptr;

}


std::unique_ptr<kgl::Variant> kgl::Variant::cloneNullVariant() const {

  StringDNA5 reference1_str(reference().getAlphabetString());
  StringDNA5 reference2_str(reference().getAlphabetString());
  VariantEvidence null_evidence; // no evidence is passed through.

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId(),
                                                                  identifier(),
                                                                  std::move(reference1_str),
                                                                  std::move(reference2_str),
                                                                  null_evidence));

  return variant_ptr;

}




// Clone with modified phase.
std::unique_ptr<kgl::Variant> kgl::Variant::clonePhase(VariantPhase phaseId) const {

  StringDNA5 reference_str(reference().getAlphabetString());
  StringDNA5 alternate_str(alternate().getAlphabetString());
  std::string variant_ident = identifier();

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId,
                                                                  std::move(variant_ident),
                                                                  std::move(reference_str),
                                                                  std::move(alternate_str),
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

  // Check using suffix and prefix
  size_t suffix_size = commonSuffix();
  size_t prefix_size = commonPrefix();
  int64_t diff_size = referenceSize() - (suffix_size + prefix_size);


    // Obvious SNP.
  if (reference_.length() == 1 and alternate_.length() == 1) {


    if (diff_size != 1) {

      ExecEnv::log().warn("Variant::isSNP; Variant: {}, Cigar: {}, prefix size: {}, suffix_size: {}, diff_size: {} not equal 1"
          , HGVS(), alternateCigar(), prefix_size, suffix_size, diff_size);

    }
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

// Check longer reference and alternate for effective SNP if a cigar of 1'X' and n'M'.
//  return ParseVCFCigar::isSNP(reference().getSequenceAsString(), alternate().getSequenceAsString());

  if (diff_size != 1) {

    ExecEnv::log().warn("Variant::isSNP; Variant: {}, Cigar: {}, prefix size: {}, suffix_size: {}, diff_size: {} not equal 1"
                        , HGVS(), alternateCigar(), prefix_size, suffix_size, diff_size);

  }

  return true;

}

// Remove common prefix and suffix nucleotides from the reference and alternate.
// .first is the trimmed reference, .second is the trimmed alternate.
std::pair<kgl::DNA5SequenceLinear, kgl::DNA5SequenceLinear> kgl::Variant::trimmedSequences() const {

  size_t prefix_size = commonPrefix();
  size_t suffix_size = commonSuffix();

  auto trimmed_reference = reference().midSequence(prefix_size, suffix_size);
  auto trimmed_alternate = alternate().midSequence(prefix_size, suffix_size);

  return { DNA5SequenceLinear(std::move(trimmed_reference)), DNA5SequenceLinear(std::move(trimmed_alternate)) };

}


// The extentOffset() of the variant is used to assess if a variant modifies a particular region of a sequence in the interval [a, b).
// With SNP variants the extentOffset().first is c = (offset() + alleleOffset()) so for "5M1X" the extent offset will be c = (offset() + 5).
// The extent size (extentOffset().second) will be 1. Thus extentOffset() will return the pair [offset()+5,1].
// The delete variant "5M4D" will have an extentOffset().first of (offset() + 5) and an extent size of 4 (reference.length() - 5).
// Thus the delete variant "5M4D" can modify the sequence [a, b) for example, (offset() + alleleOffset()) = a-2.
// The deleted nucleotides will be {a-2. a-1, a,  a+1}, the delete variant need not be contained within [a, b).
// The insert variant "5M6I", the same as an SNP, will have an extentOffset() of 5 and an extent size of 1.
// The insert variant will only modify [a, b) if (offset() + alleleOffset()) is in [a, b).
// These variants will modify a sequence [a, b) (not just translate it's offsets) if the following condition is met:
// bool modified = (extentOffset().first + extentOffset.second) > a or (extentOffset().first < b);
std::pair<kgl::ContigOffset_t, kgl::ContigSize_t> kgl::Variant::extentOffset() const {

  if (isSNP()) {

    return { offset()+alleleOffset(), 1 };

  } else {

    if (referenceSize() < alternateSize()) {
      // An insert
      return {offset()+alleleOffset(), 1 };

    } else if (referenceSize() > alternateSize()) {
      // A delete
      return { offset()+alleleOffset(), referenceSize()-alleleOffset() };

    } else {

      ExecEnv::log().warn("Variant::extentOffset; Unexpected variant: {}, cigar:{}", HGVS(), alternateCigar());

    }

  }

  return {offset()+alleleOffset(), 1 };

}


// For a sequence on the interval [a, b), Given a start offset a and a size (b-a). Determine if the variant will
// modify the sequence. Note that this is different to just translating the sequence offsets. Any upstream indel will
// modify the sequence [a, b) offsets but may not actually modify any of the nucleotides in the sequence.
bool kgl::Variant::sequenceModifier(ContigOffset_t sequence_start, ContigSize_t sequence_size) const {

  auto [extent_offset, extent_size] = extentOffset();
  return sequence_start < (extent_offset + extent_size) and extent_offset < (sequence_start + sequence_size);

}

std::string kgl::Variant::alternateCigar() const {

  return ParseVCFCigar::generateCigar(reference().getSequenceAsString(), alternate().getSequenceAsString());

}


size_t kgl::Variant::alternateSize(size_t reference_size) const {

  CigarVector cigar_vector = ParseVCFCigar::generateEditVector(reference().getSequenceAsString(), alternate().getSequenceAsString());
  return ParseVCFCigar::alternateCount(reference_size, cigar_vector);

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
