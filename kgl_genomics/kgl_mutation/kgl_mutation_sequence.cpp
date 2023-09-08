//
// Created by kellerberrin on 7/09/23.
//

#include "kgl_mutation_sequence.h"


namespace kgl = kellerberrin::genome;



bool kgl::AdjustedSequence::updateSequence() {

  for (auto const& [offset, interval_update] : interval_modify_map_) {

    switch(interval_update.variantPtr()->variantType()) {

      case VariantType::SNP:
        break;

      case VariantType::INDEL_DELETE:
        break;

      case VariantType::INDEL_INSERT:
        break;

    }

  }

  return true;

}
