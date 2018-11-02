//
// Created by kellerberrin on 26/10/18.
//

#ifndef KGL_SEQUENCE_COMPLEXITY_H
#define KGL_SEQUENCE_COMPLEXITY_H

#include "kgl_sequence_base.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class SequenceComplexity {

public:

  SequenceComplexity() = delete;
  ~SequenceComplexity() = delete;

  static double cumulativeEntropy(std::shared_ptr<const DNA5SequenceLinear> sequence, size_t word_size);

  // Calculate Shannon Entropy
  static double shannonEntropy(std::shared_ptr<const DNA5SequenceLinear> sequence);

  static size_t complexityLempelZiv(std::shared_ptr<const DNA5SequenceLinear> sequence);
  // Calculate GC content.
  static double propGC(std::shared_ptr<const DNA5SequenceLinear> sequence);


private:


};





}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_SEQUENCE_COMPLEXITY_H
