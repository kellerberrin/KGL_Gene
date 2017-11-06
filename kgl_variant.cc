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

      coding_sequences_ptr_ = findOffsetCDS(contig_offset_, gene_membership_);
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


// Returns all the coding sequences of a variant.
std:: shared_ptr<kgl::CodingSequenceArray> kgl::VariantSequence::findOffsetCDS(ContigOffset_t contig_offset,
                                                                              const GeneVector& gene_ptr_vec) {

  std:: shared_ptr<CodingSequenceArray> merged_sequences(std::make_shared<CodingSequenceArray>());

  for (const auto& gene_ptr : gene_ptr_vec) {

    std:: shared_ptr<const CodingSequenceArray> coding_sequences = kgl::GeneFeature::getCodingSequences(gene_ptr);
    merged_sequences->mergeArrays(kgl::GeneFeature::getOffsetSequences(contig_offset, coding_sequences));

  } // for each gene_ptr

  return merged_sequences;

}
