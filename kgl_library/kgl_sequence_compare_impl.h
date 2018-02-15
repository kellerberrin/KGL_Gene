//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_COMPARE_IMPL_H
#define KGL_SEQUENCE_COMPARE_IMPL_H


#include <memory>
#include <string>
#include <map>
#include "kgl_logging.h"
#include "kgl_genome_types.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class SequenceManipulation {

public:

  explicit SequenceManipulation();
  ~SequenceManipulation();

  CompareScore_t MyerHirschberg(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareDistance_t LevenshteinGlobal(const std::string& sequenceA, const std::string& sequenceB) const;

  CompareDistance_t LevenshteinLocal(const std::string& sequenceA, const std::string& sequenceB) const;

private:

  class SequenceManipImpl;       // Forward declaration of the Sequence Manipulation implementation class
  std::unique_ptr<SequenceManipImpl> sequence_manip_impl_ptr_;    // Sequence Manipulation PIMPL


};


}   // namespace genome
}   // namespace kellerberrin









#endif //KGL_KGL_SEQUENCE_COMPARE_IMPL_H
