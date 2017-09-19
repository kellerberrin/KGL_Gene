// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
//
// Created by kellerberrin on 19/09/17.
//

#ifndef KGL_LOCK_H
#define KGL_LOCK_H

#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>
#include "kgl_genome_types.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Only lock small sections of the data array for write access.
// This is a classic trade-off between the space taken by the mutex array and fine-grained speedup of access
// Advantages: Can be used on any processor and underlying count type. Disadvantage : slow.
template <typename CountType, std::size_t Granularity>
class GranularityMutex {

public:

  explicit GranularityMutex(const std::size_t array_size): array_size_(array_size), granularity_(Granularity) {

    auto lock_count = std::size_t(array_size / granularity_) + 1;
    lock_array_ = std::move(std::vector<std::mutex>{lock_count});

  }

  inline void incrementCount(CountType& counter, std::size_t array_index) {

    acquire(array_index);
    ++counter;
    release(array_index);

  }

private:

  inline void acquire(std::size_t array_index) {

    auto idx = std::size_t(array_index / granularity_);
    lock_array_[idx].lock();

  }

  inline void release(std::size_t array_index) {

    auto idx = std::size_t(array_index / granularity_);
    lock_array_[idx].unlock();

  }

  const std::size_t array_size_;
  const std::size_t granularity_;
  mutable std::vector<std::mutex> lock_array_;

};

// Does not use mutexes. Uses Intel assembler (xaddl) to implement fast locks.
// Advantage: fast. Disadvantages: Only valid on Intel processors and NucleotideReadCount_t (unsigned long) counters.
class X86Mutex {

public:

  X86Mutex(const std::size_t array_size) {}

  inline void incrementCount(NucleotideReadCount_t& counter, std::size_t array_index) {

    int inc = 1;

    __asm__ volatile("lock; xaddl %0, %1"
    : "+r" (inc), "+m" (counter) // input+output
    : // memory and condition codes changed
    : "memory", "cc"
    );

  }


};

// Use with single threaded access to the read count data..
template <typename CountType>
class NullMutex {

public:

  explicit NullMutex(const std::size_t array_size) {}

  inline void incrementCount(CountType& counter, std::size_t array_index) { ++counter; }

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_LOCK_H
