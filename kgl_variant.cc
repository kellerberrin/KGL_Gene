//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kgl_patterns.h"
#include "kgl_filter.h"
#include "kgl_variant_db.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant offset output convention.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::offsetOutput(kgl::ContigOffset_t offset, kgl::VariantOutputIndex output_base) {

  std::stringstream ss;

  switch(output_base) {

    case VariantOutputIndex::START_0_BASED:
      ss << offset;
      break;

    case VariantOutputIndex::START_1_BASED:
      ss << (offset + 1);
      break;

  }

  return ss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::VariantSequence::typestr() const {

  switch(type()) {

    case VariantSequenceType::CDS_CODING : return "Coding";

    case VariantSequenceType::INTRON : return "Intron";

    case VariantSequenceType::NON_CODING : return "NonCoding";

  }

  return "ERROR"; // Should not happen.

}

std::string kgl::VariantSequence::genomeOutput(char delimiter, VariantOutputIndex output_index) const {

  std:: stringstream ss;
// Contig.
  ss << variantSource() << delimiter <<contig()->contigId()
     << delimiter << offsetOutput(offset(), output_index)
     << delimiter << typestr() << delimiter;

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

// Convenience routine to reduceNode downstream boiler plate.
bool kgl::VariantSequence::codonOffset(ContigOffset_t& codon_offset, ContigSize_t& base_in_codon)const {

  if (not codingSequences().empty()) {

    return contig()->sequence().codonOffset(codingSequences().getFirst(), offset(), codon_offset, base_in_codon);

  } else {

    return false;

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string kgl::Variant::name() const {

 switch(variantType()) {

   case VariantType::SNP: return "SNP";

   case VariantType::DELETE: return "DEL";

   case VariantType::INSERT: return "INS";

   case VariantType::COMPOUND_SNP: return "CSNP";

   case VariantType::COMPOUND_DELETE: return "CDEL";

   case VariantType::COMPOUND_INSERT: return "CINS";

 }

 return "NOT_IMPLEMENTED";  // Not reached, to keep the compiler happy.

}


