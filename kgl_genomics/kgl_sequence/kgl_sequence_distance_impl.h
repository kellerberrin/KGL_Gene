//
// Created by kellerberrin on 22/02/18.
//

#ifndef KGL_SEQUENCE_DISTANCE_IMPL_H
#define KGL_SEQUENCE_DISTANCE_IMPL_H


#include <memory>
#include <string>
#include <map>
#include "kel_logging.h"
#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization level namespace



class SequenceDistanceImpl {

public:

  explicit SequenceDistanceImpl();
  ~SequenceDistanceImpl();

  // Distance

  [[nodiscard]] CompareDistance_t LevenshteinGlobal(const std::string& sequenceA, const std::string& sequenceB) const;

  [[nodiscard]] CompareDistance_t LevenshteinLocal(const std::string& sequenceA, const std::string& sequenceB) const;

  [[nodiscard]] CompareDistance_t globalblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;

  [[nodiscard]] CompareDistance_t localblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;

private:

  class SequenceManipImpl;       // Forward declaration of the Sequence Manipulation implementation class
  std::unique_ptr<SequenceManipImpl> sequence_manip_impl_ptr_;    // Sequence Manipulation PIMPL


};




}   // end namespace




#endif //KGL_SEQUENCE_DISTANCE_IMPL_H
