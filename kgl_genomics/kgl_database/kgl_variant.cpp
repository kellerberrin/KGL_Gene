//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kel_patterns.h"
#include "kgl_variant_factory_vcf_parse_cigar.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


[[nodiscard]] std::unique_ptr<kgl::Variant> kgl::Variant::clone() const {

  StringDNA5 reference_str(reference().getAlphabetString());
  StringDNA5 alternate_str(alternate().getAlphabetString());
  std::string variant_ident = identifier();

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId(),
                                                                  std::move(variant_ident),
                                                                  std::move(reference_str),
                                                                  std::move(alternate_str),
                                                                  evidence(),
                                                                  passFilters()));

  return variant_ptr;

}


[[nodiscard]] std::unique_ptr<kgl::Variant> kgl::Variant::cloneNullVariant() const {

  StringDNA5 reference1_str(reference().getAlphabetString());
  StringDNA5 reference2_str(reference().getAlphabetString());
  VariantEvidence null_evidence; // no evidence is passed through.

  std::unique_ptr<Variant> variant_ptr(std::make_unique<Variant>( contigId(),
                                                                  offset(),
                                                                  phaseId(),
                                                                  identifier(),
                                                                  std::move(reference1_str),
                                                                  std::move(reference2_str),
                                                                  null_evidence,
                                                                  passFilters()));

  return variant_ptr;

}



std::string kgl::Variant::typeText() const {

  switch(variantType()) {

    case VariantType::INDEL: return "INDEL";

    case VariantType::TRANSVERSION: return "TRANSVERSION";

    case VariantType::TRANSITION: return "TRANSITION";

  }

  return "NOT_IMPLEMENTED";  // Not reached, to keep the compiler happy.

}


kgl::VariantType kgl::Variant::variantType() const {

  if (not isSNP()) {

    return VariantType::INDEL;

  } else {

    if (DNA5::isTransition(reference().at(0), alternate().at(0))) {

      return VariantType::TRANSITION;

    } else {

     return VariantType::TRANSVERSION;

    }

  }

}


std::string kgl::Variant::genomeOutput(char delimiter, VariantOutputIndex output_index) const {

  std:: stringstream ss;
  // Contig.
  ss << contigId() << delimiter;
  if (phaseId() == VariantPhase::UNPHASED) {

    ss << "Unphased" << delimiter;

  } else {

    ss << "Phase:" << static_cast<size_t>(phaseId()) << delimiter;

  }
  ss << offsetOutput(offset(), output_index) << delimiter;

  return ss.str();

}



std::string kgl::Variant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << identifier() << delimiter << typeText() << delimiter << alternateSize() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence().output(delimiter, output_index);

  }

  return ss.str();

}

// Identifier is not used in these equivalence relationships.
// Equal Phase
bool kgl::Variant::equivalent(const Variant& cmp_var) const {

  return contigId() == cmp_var.contigId()
         and phaseId() == cmp_var.phaseId()
         and offset() == cmp_var.offset()
         and variantType() == cmp_var.variantType()
         and reference().getSequenceAsString() == cmp_var.reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_var.alternate().getSequenceAsString();

}

// Opposite Phase
bool kgl::Variant::homozygous(const Variant& cmp_var) const {

  return contigId() == cmp_var.contigId()
         and phaseId() != cmp_var.phaseId()
         and offset() == cmp_var.offset()
         and variantType() == cmp_var.variantType()
         and reference().getSequenceAsString() == cmp_var.reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_var.alternate().getSequenceAsString();

}

// No Phase
bool kgl::Variant::analogous(const Variant& cmp_var) const {

  return contigId() == cmp_var.contigId()
         and offset() == cmp_var.offset()
         and variantType() == cmp_var.variantType()
         and reference().getSequenceAsString() == cmp_var.reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_var.alternate().getSequenceAsString();

}




// Order variant types.
bool kgl::Variant::lessThan(const Variant& cmp_var) const {


  if (contigId() < cmp_var.contigId()) {

    return true;

  } else if (contigId() > cmp_var.contigId()) {

    return false;

  } else if (phaseId() < cmp_var.phaseId()) {

    return true;

  } else if (phaseId() > cmp_var.phaseId()) {

    return false;

  } else if (offset() < cmp_var.offset()) {

    return true;

  } else if (offset() > cmp_var.offset()) {

    return false;

  }

  if (reference().getSequenceAsString() < cmp_var.reference().getSequenceAsString()) {

    return true;

  } else if (reference().getSequenceAsString() > cmp_var.reference().getSequenceAsString()) {

    return false;

  } else if (alternate().getSequenceAsString() < cmp_var.alternate().getSequenceAsString()) {

    return true;

  }

  return false;

}


std::string kgl::Variant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << reference().getSequenceAsString() << ">" << offsetOutput(offset(), output_index) << ">";
  ss << alternate().getSequenceAsString() << delimiter;
  ss << alternateCigar() << delimiter;

  return ss.str();

}



std::string kgl::Variant::alternateCigar() const {

  return ParseVCFCigar::generateCigar(reference().getSequenceAsString(), alternate().getSequenceAsString());

}


size_t kgl::Variant::alternateSize(size_t reference_size) const {

  CigarVector cigar_vector = ParseVCFCigar::generateEditVector(reference().getSequenceAsString(), alternate().getSequenceAsString());
  return ParseVCFCigar::alternateCount(reference_size, cigar_vector);

}



// Unique up to phase
std::string kgl::Variant::alleleHash() const {

  return reference().getSequenceAsString() + ":" + alternate().getSequenceAsString();

}

// Phase specific hash
std::string kgl::Variant::allelePhaseHash() const {

  return std::to_string(static_cast<uint8_t>(phaseId())) + ":" + reference().getSequenceAsString() + ":" + alternate().getSequenceAsString();

}


// Unique upto phase.
std::string kgl::Variant::variantHash() const {

  return contigId() + ":"  + std::to_string(offset()) + ":" + alleleHash();

}

// Phase specific hash
std::string kgl::Variant::variantPhaseHash() const {

  return contigId() + ":"  + std::to_string(offset()) + ":" + allelePhaseHash();

}
