//
// Created by kellerberrin on 15/05/18.
//

#include "kgd_combinatorial.h"

namespace kgd = kellerberrin::deconvolv;


void kgd::BinaryPermutations::checkBitSize() {

  if (bits_ > bitmask_.size()) {

    ExecEnv::log().critical("BinaryPermutations() max bits: {}, requested: {} bits", bitmask_.size(), bits_);

  }


  bitmask_ = std::numeric_limits<std::size_t>::max();

  size_t shift_right = bitmask_.size() - bits_;

  bitmask_ >>= shift_right;

  max_index_ = bitmask_.to_ullong();

}

