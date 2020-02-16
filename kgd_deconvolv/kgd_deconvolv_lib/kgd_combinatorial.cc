//
// Created by kellerberrin on 15/05/18.
//

#include "kgd_combinatorial.h"

namespace kgd = kellerberrin::deconvolv;




std::vector<std::vector<size_t> > kgd::BinaryPermutations::enumerateBinaryMatrix() const {
  // This function enumerate all possible binary combinations of k elements
  //TODO: Exponential explosion causes memory exhaustion for (k > 5). Replace/remove this logic.

  std::vector<std::vector<size_t> > ret;

  for (size_t n = 0; n < bitArraySize(); ++n) {

    ret.push_back(reverseBitArray(n));

  }

  return ret;

}


std::vector<size_t> kgd::BinaryPermutations::reverseBitArray(size_t n) const {

  std::vector<size_t> ret;

  while (n > 0) {

    size_t binary_digit = (n & 0x1u) ? 1 : 0;

    ret.push_back(binary_digit);

    n = n >> 0x1u;

  }

  while(ret.size() < bits_) {

    ret.push_back(0);

  }

  if (ret.size() > bits_) {

    ExecEnv::log().critical("BinaryPermutations; Vector Binary length: {} greater than specified: {}", ret.size(), bits_);

  } else if (ret.size() < bits_){

    ExecEnv::log().warn("BinaryPermutations; Vector Binary length: {} less than specified: {}", ret.size(), bits_);

  }

  std::reverse(ret.begin(), ret.end()); // reverse vector e.g. {1, 0, 1, 1, 0} (20) to { 0, 1, 1, 0, 1 } (13)

  return ret;

}


