//
// Created by kellerberrin on 29/10/17.
//


#ifndef KGL_MUTANT_SEQUENCE_H
#define KGL_MUTANT_SEQUENCE_H

#include "kgl_base_sequence.h"
#include "kgl_variant.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class MutantSequence   {

public:

  explicit MutantSequence(std::shared_ptr<DNA5Sequence> sequence_ptr) : sequence_ptr_(sequence_ptr) {};
  MutantSequence() = delete;
  virtual ~MutantSequence() = default;


  std::shared_ptr<DNA5Sequence> mutantSequence(std::shared_ptr<const DNA5Sequence> sequence_ptr,
                                                  const OffsetVariantMap& variant_map,
                                                  const SortedCDS& sorted_cds);


private:

  std::shared_ptr<DNA5Sequence> sequence_ptr_;

};



}   // namespace genome
}   // namespace kellerberrin






#endif //KGL_MUTANT_SEQUENCE_H
