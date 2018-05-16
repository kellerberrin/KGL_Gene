//
// Created by kellerberrin on 15/05/18.
//

#ifndef KGD_COMBINATORIAL_H
#define KGD_COMBINATORIAL_H

#include <vector>
#include <bitset>
#include <limits>
#include "kgd_deconvolv_app.h"

namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


/// Max 64 bits.
using BitVector = std::bitset<std::numeric_limits<size_t>::digits>;

/// This class removes the need to generate a vector of binary bit patterns by
/// returning the bit pattern for any arbitrary index.
class BinaryPermutations {

public:

  BinaryPermutations(size_t bits) : bits_(bits) { checkBitSize(); }
  ~BinaryPermutations() = default;

  size_t bitSize() const { return bits_; }

  size_t maxIndex() const { return max_index_; }

  template<typename T> std::vector<T> getBitpattern(size_t index) const {

    BitVector bit_index = index;

    if (index > maxIndex()) {

      bit_index &= bitmask_;
      ExecEnv::log().error("BinaryPermutations() index: {} out of bounds, masked to: {}", index, bit_index.to_ullong());

    }

    std::vector<T> bit_array;

    for (size_t idx = 0; idx < bits_; ++idx) {

      T bit_value = bit_index.test(idx) ? 1 : 0;

      bit_array.push_back(bit_value);

    }

    return bit_array;

  }

private:

  size_t bits_;
  size_t max_index_;
  BitVector bitmask_;

  void checkBitSize();

};




}   // organization level namespace
}   // project level namespace



#endif //KGL_KGD_COMBINATORIAL_H
