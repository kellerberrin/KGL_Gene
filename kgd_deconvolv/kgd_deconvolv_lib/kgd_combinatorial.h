//
// Created by kellerberrin on 15/05/18.
//

#ifndef KGD_COMBINATORIAL_H
#define KGD_COMBINATORIAL_H

#include <vector>
#include <bitset>
#include <limits>
#include "kgd_deconvolv_app.h"

namespace kellerberrin::deconvolv {    // organization level namespace



/// This class removes the need to generate a vector of binary bit patterns by
/// returning the bit pattern for any arbitrary index.
class BinaryPermutations {

public:

  explicit BinaryPermutations(size_t bits) : bits_(bits) { bit_array_size_ = 0x1u << bits_; }
  ~BinaryPermutations() = default;

  // Return a matrix of all reversed array of bits.
  std::vector<std::vector<size_t> > enumerateBinaryMatrix() const;
  // The maximum size of the bit array
  size_t bitArraySize() const { return bit_array_size_; }
  // A reverse bit array of bits length for the specified number n.
  std::vector<size_t> reverseBitArray(size_t n) const;

private:

  size_t bits_;
  size_t bit_array_size_;

};




}   // end namespace



#endif //KGL_KGD_COMBINATORIAL_H
