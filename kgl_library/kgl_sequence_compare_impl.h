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

  CompareScore_t MyerHirschbergGlobal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareScore_t MyerHirschbergLocal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareDistance_t LevenshteinGlobal(const std::string& sequenceA, const std::string& sequenceB) const;

  CompareDistance_t LevenshteinLocal(const std::string& sequenceA, const std::string& sequenceB) const;

  CompareDistance_t globalblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;

  std::string editItems(const std::string& reference, const std::string& mutant, char delimiter, VariantOutputIndex index_offset);

private:

  class SequenceManipImpl;       // Forward declaration of the Sequence Manipulation implementation class
  std::unique_ptr<SequenceManipImpl> sequence_manip_impl_ptr_;    // Sequence Manipulation PIMPL


};


}   // namespace genome
}   // namespace kellerberrin









#endif //KGL_KGL_SEQUENCE_COMPARE_IMPL_H
