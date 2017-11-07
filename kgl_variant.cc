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

  switch(genomeType()) {

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

kgl::VariantSequenceType kgl::VariantSequence::genomeType() {

  if (variant_genome_type_ == VariantSequenceType::UNKNOWN) {

    if (contig_ptr_->findGenes(contig_offset_, gene_membership_)) {

      coding_sequences_ptr_ = kgl::GeneFeature::findOffsetCDS(contig_offset_, gene_membership_);
      if(coding_sequences_ptr_->size() > 0) {

        variant_genome_type_ = VariantSequenceType::CDS_CODING;

      }
      else {

        variant_genome_type_ = VariantSequenceType::INTRON;

      }

    } else {

      variant_genome_type_ = VariantSequenceType::NON_CODING;

    }

  }

  return variant_genome_type_;

}

