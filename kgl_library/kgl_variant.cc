//
// Created by kellerberrin on 13/10/17.
//

#include <ostream>
#include "kgl_variant.h"
#include "kgl_patterns.h"
#include "kgl_filter.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Genome information of the variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


const kgl::PhaseId_t kgl::VariantSequence::UNPHASED;


std::string kgl::VariantSequence::genomeOutput(char delimiter, VariantOutputIndex output_index) const {

  std:: stringstream ss;
// Contig.
  ss << genomeId() << delimiter
     << contigId() << delimiter;
  if (phaseId() == UNPHASED) {

    ss << "UNPHASED" << delimiter;

  } else {

    ss << "PHASE:" << static_cast<size_t>(phaseId()) << delimiter;

  }
  ss << offsetOutput(offset(), output_index) << delimiter;

  return ss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Variant class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string kgl::Variant::name() const {

 switch(variantType()) {

   case VariantType::SNP: return "SNP";

   case VariantType::DELETE: return "DEL";

   case VariantType::INSERT: return "INS";

   case VariantType::VCF_VARIANT: return "VCF";

   case VariantType::COMPOUND_DELETE: return "CDEL";

   case VariantType::COMPOUND_INSERT: return "CINS";

 }

 return "NOT_IMPLEMENTED";  // Not reached, to keep the compiler happy.

}


// Relevant for compound variants.
// Used to see if variants are mutating the same section of DNA.
bool kgl::Variant::offsetOverlap(const Variant& cmp_var) const {

  return false; // test if overlapping variants generate an error.

  // First check if the variants are on the same contig.
  if (contigId() != cmp_var.contigId()) {

    return false;

  }

  // On the same contig so check for overlap.
  ContigOffset_t start = offset();
  ContigOffset_t end = offset() + size() - 1;

  ContigOffset_t cmp_start = cmp_var.offset();
  ContigOffset_t cmp_end = cmp_var.offset() + size() - 1;

  bool overlap_offset = ((start <= cmp_start) and (cmp_start <= end)) or ((start <= cmp_end) and (cmp_end <= end));

  if (overlap_offset) {

    ExecEnv::log().info("Variant::offsetOverlap(), variants overlap\n{}\n{}",
                        output(' ', VariantOutputIndex::START_0_BASED, true),
                        cmp_var.output(' ', VariantOutputIndex::START_0_BASED, true));

  }

  return overlap_offset;

}