//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kgl_patterns.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::VariantSequence::typestr() const {

  switch(type()) {

    case VariantSequenceType::UNKNOWN: return "UNKNOWN";

    case VariantSequenceType::CDS_CODING : return "CDS Coding";

    case VariantSequenceType::INTRON : return "INTRON";

    case VariantSequenceType::NON_CODING : return "NON Coding";

  }

  return "ERROR"; // Should not happen.

}


std::string kgl::VariantSequence::genomeOutput() const {

  std:: stringstream ss;
// Contig.

  ss << contig()->contigId();
  ss << " " << typestr() << " ";


  return ss.str();

}


void kgl::VariantSequence::defineVariantType(std::shared_ptr<const GeneFeature> gene_ptr,
                                             std::shared_ptr<const CodingSequence> coding_sequence_ptr)
{

  gene_membership_ = gene_ptr;
  coding_sequence_ptr_ = coding_sequence_ptr;

  if (gene_ptr != nullptr and coding_sequence_ptr != nullptr) {

    variant_genome_type_ = VariantSequenceType::CDS_CODING;

  } else if (gene_ptr != nullptr) {


    variant_genome_type_ = VariantSequenceType::INTRON;


  } else {

    variant_genome_type_ = VariantSequenceType::NON_CODING;

  }

}


