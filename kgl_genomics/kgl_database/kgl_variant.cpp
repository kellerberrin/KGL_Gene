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




std::string kgl::Variant::typeText() const {

  switch(variantType()) {

    case VariantType::INDEL: return "INDEL";

    case VariantType::TRANSVERSION: return "TRANSVERSION";

    case VariantType::TRANSITION: return "TRANSITION";

  }

  return "NOT_IMPLEMENTED";  // Not reached, to keep the compiler happy.

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
  ss << offsetOutput(offset(), output_index)
     << "(+" << alleleOffset() << ")" << delimiter ;

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


std::string kgl::Variant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << reference().getSequenceAsString() << ">" << offsetOutput(referenceOffset(), output_index) << ">";
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


// Unique upto phase.
std::string kgl::Variant::HGVS() const {

  // 1 is added to the offset to make it 1-based as is the standard in Gffs etc.
  return contigId() + ":g."  + std::to_string(referenceOffset() + 1) + reference().getSequenceAsString() + ">" + alternate().getSequenceAsString();

}

// Phase specific hash
std::string kgl::Variant::HGVS_Phase() const {

  // 1 is added to the offset to make it 1-based as is the standard in Gffs etc.
  return contigId() + ":g."  + std::to_string(referenceOffset() + 1) + reference().getSequenceAsString() +
         ">" + alternate().getSequenceAsString() + ":" + std::to_string(static_cast<uint8_t>(phaseId()));

}
