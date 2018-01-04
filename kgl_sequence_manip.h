//
// Created by kellerberrin on 4/01/18.
//

#ifndef KGL_SEQUENCE_MANIP_H
#define KGL_SEQUENCE_MANIP_H


#include <memory>
#include <string>
#include <map>
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class SequenceManipulation {

public:

  explicit SequenceManipulation();
  ~SequenceManipulation();

  std::string compareSequences(const std::string& reference_str, const std::string& compare_str) const;


private:

  class SequenceManipImpl;       // Forward declaration of the Sequence Manipulation implementation class
  std::unique_ptr<SequenceManipImpl> sequence_manip_impl_ptr_;    // Sequence Manipulation PIMPL


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_MANIP_H
