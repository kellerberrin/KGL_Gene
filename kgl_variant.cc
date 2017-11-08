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

kgl::VariantSequenceType kgl::VariantSequence::type() const {

  if (not codingSequences().empty()) {

    return VariantSequenceType::CDS_CODING;

  } else if (not geneMembership().empty()) {

    return VariantSequenceType::INTRON;

  } else {

    return VariantSequenceType::NON_CODING;

  }

}

void kgl::VariantSequence::defineIntron(std::shared_ptr<const GeneFeature> gene_ptr)
{

  if (gene_ptr) {

    coding_sequences_.getMap().clear();
    gene_membership_.clear();
    gene_membership_.push_back(gene_ptr);

  } else {


    ExecEnv::log().error("Variant contig: {} offset: {}; Attempted to define intron with null pointer",
                         contig()->contigId(), offset());

  }

}

void kgl::VariantSequence::defineCoding(std::shared_ptr<const CodingSequence> coding_sequence_ptr)
{

  if (coding_sequence_ptr) {

    coding_sequences_.getMap().clear();
    coding_sequences_.insertCodingSequence(coding_sequence_ptr);
    gene_membership_.clear();
    gene_membership_.push_back(coding_sequence_ptr->getGene());

  } else {


    ExecEnv::log().error("Variant contig: {} offset: {}; Attempted to define coding sequence with null pointer",
                         contig()->contigId(), offset());

  }

}

void kgl::VariantSequence::defineNonCoding()
{

  coding_sequences_.getMap().clear();
  gene_membership_.clear();

}



