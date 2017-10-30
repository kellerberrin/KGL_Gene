//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SEQUENCE_MUTANT_H
#define KGL_SEQUENCE_MUTANT_H


#include "kgl_sequence_base.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class MutantSequence   {

public:

  MutantSequence() = default;
  ~MutantSequence() = default;


  std::shared_ptr<DNA5Sequence> mutantSequence(std::shared_ptr<const DNA5Sequence> sequence_ptr,
                                               const OffsetVariantMap& variant_map,
                                               const SortedCDS& sorted_cds);

private:


};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SEQUENCE_MUTANT_H
