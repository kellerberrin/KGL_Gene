//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_COMPARE_H
#define KGL_SEQUENCE_COMPARE_H

#include "kgl_genome_types.h"
#include "kgl_sequence_virtual.h"
#include "kgl_sequence_compare_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class SequenceComparison {

public:

  SequenceComparison() = default;
  virtual ~SequenceComparison() = default;

  virtual std::string comparisonType() const = 0;

  virtual CompareScore_t comparison(std::shared_ptr<const VirtualSequence> sequenceA,
                          std::shared_ptr<const VirtualSequence> sequenceB,
                          std::string& comparison_string) const = 0;

private:


};


class MyerHirschbergComparison: public SequenceComparison {

public:

  MyerHirschbergComparison() = default;
  virtual ~MyerHirschbergComparison() = default;

  std::string comparisonType() const override { return "MyerHirschberg"; }


  virtual CompareDistance_t distance(std::shared_ptr<const VirtualSequence> sequenceA,
                                     std::shared_ptr<const VirtualSequence> sequenceB,
                                     CompareDistance_t& distance,
                                     std::string& comparison_string) const override {

    return SequenceManipulation().MyerHirschberg(sequenceA->getSequenceAsString(),
                                                 sequenceB->getSequenceAsString(),
                                                 comparison_string);

  }


private:


};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_COMPARE_H
