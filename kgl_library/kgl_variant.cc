//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kel_patterns.h"
#include "kgl_filter.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"
#include "kgl_variant_factory_vcf_parse_cigar.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::string kgl::VariantSequence::genomeOutput(char delimiter, VariantOutputIndex output_index) const {

  std:: stringstream ss;
// Contig.
  ss << contigId() << delimiter;
  if (phaseId() == UNPHASED) {

    ss << "Unphased" << delimiter;

  } else {

    ss << "Phase:" << static_cast<size_t>(phaseId()) << delimiter;

  }
  ss << offsetOutput(offset(), output_index) << delimiter;

  return ss.str();

}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::Variant::name() const {

  switch(variantType()) {

    case VariantType::VCF_VARIANT: return "VCF";

  }

  return "NOT_IMPLEMENTED";  // Not reached, to keep the compiler happy.

}

std::string kgl::Variant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << alternateSize() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence().output(delimiter, output_index);

  }

  return ss.str();

}



bool kgl::Variant::equivalent(const Variant& cmp_var) const {

  return contigId() == cmp_var.contigId()
         and phaseId() == cmp_var.phaseId()
         and offset() == cmp_var.offset()
         and variantType() == cmp_var.variantType()
         and reference().getSequenceAsString() == cmp_var.reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_var.alternate().getSequenceAsString();

}


bool kgl::Variant::homozygous(const Variant& cmp_var) const {

  return contigId() == cmp_var.contigId()
         and phaseId() != cmp_var.phaseId()
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

  } else if (variantType() < cmp_var.variantType()) {

    return true;

  } else if (variantType() > cmp_var.variantType()) {

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
