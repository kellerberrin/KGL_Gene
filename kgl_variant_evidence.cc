//
// Created by kellerberrin on 4/12/17.
//



#include "kgl_variant_evidence.h"


namespace kgl = kellerberrin::genome;




std::string kgl::ReadCountEvidence::output(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << mutantCount() << "/" << readCount() << delimiter;
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << ExtendDNA5::convertToChar(ExtendDNA5::offsetToNucleotide(idx)) << ":" << countArray()[idx] << delimiter;
  }

  return ss.str();

}

