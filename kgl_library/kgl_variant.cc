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

