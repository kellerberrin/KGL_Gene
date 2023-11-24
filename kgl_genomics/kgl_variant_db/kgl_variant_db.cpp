//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant_db.h"
#include "kgl_variant_filter_type.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


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

  // Check to be sure.
  if (not variant_ptr->isCanonical()) {

    ExecEnv::log().error("Variant::cloneCanonical; variant: {} is NOT canonical", variant_ptr->HGVS());

  }

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

    return (referenceSize() < alternateSize()) ?  VariantType::INDEL_INSERT : VariantType::INDEL_DELETE;

  } else {

    return VariantType::SNP;

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
// Canonical SNPs are represented as '1X', deletes as '1MnD' and inserts as '1MnI'.
// The first argument is the adjusted reference.
// The second is the adjusted alternate.
// The third is the adjusted variant offset.
// The adjusted variant offset always indicates the offset
// at which the first nucleotide of the canonical reference and canonical alternate is found.
std::tuple<kgl::DNA5SequenceLinear, kgl::DNA5SequenceLinear, kgl::ContigOffset_t> kgl::Variant::canonicalSequences() const {

  // Variant may already be in canonical form.
  if (isCanonical()) {

    return { DNA5SequenceLinear(StringDNA5(reference().getAlphabetString())),
             DNA5SequenceLinear(StringDNA5(alternate().getAlphabetString())),
             offset() };

  }

  // Else modify prefix and suffix subsequences to canonical form.
  size_t prefix_size = commonPrefix();
  prefix_size = prefix_size > 0 ? (prefix_size - 1) : 0;  // Adjust for '1MnD' and '1MnI'.
  size_t suffix_size = commonSuffix();
  size_t min_size = std::min(referenceSize(), alternateSize());
  int64_t adj_suffix_size = std::min(min_size - prefix_size - 1, suffix_size); // Adjust for '1MnD' and '1MnI'.
  adj_suffix_size = adj_suffix_size < 0 ? 0 : adj_suffix_size;

  auto canonical_reference = reference().removePrefixSuffix(prefix_size, adj_suffix_size);
  auto canonical_alternate = alternate().removePrefixSuffix(prefix_size, adj_suffix_size);
  ContigOffset_t canonical_offset = offset() + prefix_size;

  return { DNA5SequenceLinear(std::move(canonical_reference)),
           DNA5SequenceLinear(std::move(canonical_alternate)),
           canonical_offset };

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

// The interval() of the variant is used to assess if a CANONICAL variant modifies a particular interval [a, b).
// An SNP offset is an interval size 1 with lower() = offset().
// A delete interval is the number of deleted nucleotides with lower() = (offset() + 1).
// An insert interval is the number of inserted nucleotides with lower= (offset() + 1).
std::pair<kgl::VariantType, kel::OpenRightUnsigned> kgl::Variant::modifyInterval() const {

  auto const [canonical_ref, canonical_alt, canonical_offset] = canonicalSequences();

  switch(variantType()) {

    case VariantType::INDEL_INSERT: {

      size_t insert_size = canonical_alt.length() - canonical_ref.length();
      ContigOffset_t insert_offset = canonical_offset + canonical_ref.length();
      OpenRightUnsigned insert_interval{insert_offset, insert_offset + insert_size};
      return {VariantType::INDEL_INSERT, insert_interval};

    }

    case VariantType::INDEL_DELETE: {

      size_t delete_size = canonical_ref.length() - canonical_alt.length();
      ContigOffset_t delete_offset = canonical_offset + canonical_alt.length();
      OpenRightUnsigned delete_interval {delete_offset, delete_offset + delete_size};
      return {VariantType::INDEL_DELETE, delete_interval};

    }

    default:
    case VariantType::SNP: {

      OpenRightUnsigned snp_interval{canonical_offset, canonical_offset + canonical_ref.length()};
      return {VariantType::SNP, snp_interval};

    }

  }

}

// Modify the insert interval to [member_interval.lower()-1, member_interval.lower+1)
std::pair<kgl::VariantType, kel::OpenRightUnsigned> kgl::Variant::memberInterval() const {

  auto [variant_type, member_interval] = modifyInterval();
  // An insert variant must modify the specified region.
  // Therefore the INSERT membership interval must be contained in the modified interval.
  // Otherwise the INSERTed interval is adjacent to the modified interval and does not modify it.
  if (variant_type == VariantType::INDEL_INSERT) {

    member_interval.resize(member_interval.lower()-1, member_interval.lower()+1);

  }

  return {variant_type, member_interval};

}


std::string kgl::Variant::cigar() const {

  return ParseVCFCigar::generateCigar(std::string(reference().getStringView()), std::string(alternate().getStringView()));

}

// Unique upto phase.
std::string kgl::Variant::HGVS() const {

  return std::format("{}:g.{}{}>{}", contigId(), offset(), reference().getStringView(), alternate().getStringView());

}

// Phase specific hash
std::string kgl::Variant::HGVS_Phase() const {

  return std::format("{}:g.{}{}>{}:{}", contigId(), offset(), reference().getStringView(), alternate().getStringView(), static_cast<uint8_t>(phaseId()));

}
