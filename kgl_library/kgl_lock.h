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

class AutoMutex {

public:

  AutoMutex(std::mutex& mutex) : mutex_(mutex) { mutex_.lock(); }
  AutoMutex(const AutoMutex&) = delete;
  ~AutoMutex() { mutex_.unlock(); }

  const AutoMutex& operator=(const AutoMutex&) = delete;

private:

  std::mutex& mutex_;

};


// Provides a mutex for every element of an array. Classic speed/space trade off.
// Advantages: Can be used on any processor and underlying data type. Disadvantages : slow and memory intensive.
// Is significantly faster, and uses significantly more memory, than the GranularityLock defined below .
class ArrayMutex {

public:

  explicit ArrayMutex(const std::size_t array_size): array_size_(array_size),
                                                     lock_array_ptr_(std::make_unique<std::mutex[]>(array_size)) {}


  inline void acquire(std::size_t array_index) {

    lock_array_ptr_[array_index].lock();

  }

  inline void release(std::size_t array_index) {

    lock_array_ptr_[array_index].unlock();

  }

private:

  const std::size_t array_size_;
  mutable std::unique_ptr<std::mutex[]> lock_array_ptr_;

};


// Only lock small sections of the array for write access.
// This is a classic trade-off between the space taken by the mutex array and fine-grained speedup of access
// Advantages: Can be used on any processor and underlying data type. Disadvantage : slow.
// Will use significantly less memory and be slower than the ArrayMutex defined above.
template <std::size_t Granularity = 1000>
class GranularityMutex {

public:

  explicit GranularityMutex(const std::size_t array_size): array_size_(array_size), granularity_(Granularity) {

    auto lock_count = std::size_t(array_size / granularity_) + 1;
    lock_array_ptr_.reset(new std::mutex[lock_count]);

  }

  inline void acquire(std::size_t array_index) {

    auto idx = std::size_t(array_index / granularity_);
    lock_array_ptr_[idx].lock();

  }

  inline void release(std::size_t array_index) {

    auto idx = std::size_t(array_index / granularity_);
    lock_array_ptr_[idx].unlock();

  }

private:

  const std::size_t array_size_;
  const std::size_t granularity_;
  mutable std::unique_ptr<std::mutex[]> lock_array_ptr_;

};

// Only for single thread access.
// Never use this as a multi-thread strategy when updating inserted sequences.
// A segmentation fault and sorrow will surely follow.
class NullMutex {

public:

  explicit NullMutex(const std::size_t) {}

  inline void acquire(std::size_t) {}

  inline void release(std::size_t) {}

};


// Nucleotide Count Locking strategies.
// If the nucleotide count array is not locked then some counts will be lost.

// Use a mutex strategy, potentially slow and memory intensive (ArrayMutex) but very general.
template <class LockStrategy = ArrayMutex, typename CountType = NucleotideReadCount_t>
class MutexCountLock : protected LockStrategy {

public:

  explicit MutexCountLock(const std::size_t array_size): LockStrategy(array_size) {}

  inline void incrementCount(CountType& counter, std::size_t array_index) {

    LockStrategy::acquire(array_index);
    ++counter;
    LockStrategy::release(array_index);

  }


};

// Does not use mutexes. Uses Intel assembler (xaddl) to implement fast locks.
// Advantage: fast. Disadvantages: Only valid on Intel processors and NucleotideReadCount_t (unsigned long) counters.
class X86CountLock {

public:

  explicit X86CountLock(const std::size_t) {}

  inline void incrementCount(NucleotideReadCount_t& counter, std::size_t) {

    int inc = 1;

    __asm__ volatile("lock; xaddl %0, %1"
    : "+r" (inc), "+m" (counter) // input+output
    : // memory and condition codes changed
    : "memory", "cc"
    );

  }

};

// Only use with single thread access.
// With multi-thread nucleotide count updates, this will be fast, but you will lose some (small number) counts.
template <typename CountType = NucleotideReadCount_t>
class NullCountLock {

public:

  explicit NullCountLock(const std::size_t array_size) {}

  inline void incrementCount(CountType& counter, std::size_t array_index) { ++counter; }

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_LOCK_H
