//
// Created by kellerberrin on 31/08/18.
//


#include "kgl_variant_vcf.h"
#include "kgl_sequence_offset.h"
#include "kgl_variant_factory_vcf_parse_cigar.h"

#include <algorithm>


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::string kgl::VCFVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence().output(delimiter, output_index);

  }

  return ss.str();

}



bool kgl::VCFVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const VCFVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return genomeId() == cmp_snp->genomeId()
         and contigId() == cmp_snp->contigId()
         and phaseId() == cmp_snp->phaseId()
         and offset() == cmp_snp->offset()
         and variantType() == cmp_snp->variantType()
         and reference().getSequenceAsString() == cmp_snp->reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_snp->alternate().getSequenceAsString();

}


// Order variant types.
bool kgl::VCFVariant::lessThan(const Variant& cmp_var) const {


  if (genomeId() < cmp_var.genomeId()) {

    return true;

  } else if (contigId() < cmp_var.contigId()) {

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

  auto cmp_snp = dynamic_cast<const VCFVariant*>(&cmp_var);

  if (not cmp_snp) {

    // Must be a variant type == snp type.
    ExecEnv::log().error("SNPVariant::lessThan; Expected VCFVariant, got: {}", cmp_var.output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  if (reference().getSequenceAsString() < cmp_snp->reference().getSequenceAsString()) {

    return true;

  } else if (reference().getSequenceAsString() > cmp_snp->reference().getSequenceAsString()) {

    return false;

  } else if (alternate().getSequenceAsString() < cmp_snp->alternate().getSequenceAsString()) {

    return true;

  }

  return false;

}


std::string kgl::VCFVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << reference().getSequenceAsString() << ">" << offsetOutput(offset(), output_index) << ">";
  ss << alternate().getSequenceAsString() << delimiter;
  ss << alternateCigar() << delimiter;

  return ss.str();

}



std::string kgl::VCFVariant::alternateCigar() const {

  return ParseVCFCigar::generateCigar(reference().getSequenceAsString(), alternate().getSequenceAsString());

}


size_t kgl::VCFVariant::alternateSize(size_t reference_size) const {

  CigarVector cigar_vector = ParseVCFCigar::generateEditVector(reference().getSequenceAsString(), alternate().getSequenceAsString());
  return ParseVCFCigar::alternateCount(reference_size, cigar_vector);

}
