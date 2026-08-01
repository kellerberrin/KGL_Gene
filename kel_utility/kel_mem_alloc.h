//
// Created by kellerberrin on 14/2/21.
//

#ifndef KEL_MEM_ALLOC_H
#define KEL_MEM_ALLOC_H


#include <memory_resource>
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstring>
#include <malloc.h>
#include <new>

namespace kellerberrin {   //  organization level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A heap memory auditing object.
// These functions are thread safe in linux.
// Because block size is added to the front of these memory allocations newMem(), newArray()
// must be paired with deleteMem(), deleteArray() respectively.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AuditMemory  {

public:

  AuditMemory() = delete;
  ~AuditMemory() = delete;

  // Force free store release using malloc trim.
  // Malloc trim is a GNU/glibc-specific extension.
  // It will not compile on macOS, Windows, and some non-glibc Linux systems.
  static void trimFreeStore() {
#ifdef __GLIBC__
    malloc_trim(0);
#endif
  }

  // Audited alloc functions.
  template<class T>
  [[nodiscard]] static T* newMem(std::size_t mem_size);
  template<class T>
  [[nodiscard]] static T* newArray(std::size_t array_size) { return newMem<T>(sizeof(T) * array_size); }
  // Audited free functions.
  template<class T>
  static void deleteMem(T* mem_ptr);
  template<class T>
  static void deleteArray(T* mem_ptr) { deleteMem(mem_ptr); }

  template<class T>
  [[nodiscard]] static size_t alignedArray(size_t array_size) { return alignedSize(sizeof(T) *  array_size); }

  [[nodiscard]] static size_t deallocatedBytes() { return deallocated_bytes_; }
  [[nodiscard]] static size_t allocatedBytes() { return allocated_bytes_; }
  [[nodiscard]] static size_t allocations() { return allocations_; }
  [[nodiscard]] static size_t deallocations() { return deallocations_; }
  [[nodiscard]] static size_t countMaxAlign() { return count_max_align_; }
  [[nodiscard]] static size_t additionalAlignBytes() { return additional_align_bytes_; }
  [[nodiscard]] static size_t alignedSize(size_t mem_size);

private:

  // Threadsafe malloc and free statistics to aid in the search for memory leaks.
  inline static std::atomic<size_t> deallocated_bytes_{0};
  inline static std::atomic<size_t> allocated_bytes_{0};
  inline static std::atomic<size_t> allocations_{0};
  inline static std::atomic<size_t> deallocations_{0};
  inline static std::atomic<size_t> count_max_align_{0};
  inline static std::atomic<size_t> additional_align_bytes_{0};

   // Number of bytes to allocate for a memory aligned size_t.
  constexpr static size_t SIZE_ALIGN_BYTES = std::max(alignof(std::max_align_t), sizeof(std::size_t));

};

// Warning Will Robinson! Nasty pointer munging ahead.
template<class T>
T* AuditMemory::newMem(std::size_t mem_size)
{

  if (mem_size == 0) {

    ++mem_size; // avoid std::malloc(0) which may return nullptr on success

  }

  size_t aligned_size = alignedSize(mem_size);;
  if (aligned_size == mem_size) {

    ++count_max_align_;

  } else {

    additional_align_bytes_ += aligned_size - mem_size;
    mem_size = aligned_size;

  }

  // Thread safe in linux.
  void *ptr = malloc(SIZE_ALIGN_BYTES + mem_size);
  if (ptr != nullptr) {

    // Don't count ALIGN_SIZE
    allocated_bytes_ += mem_size;
    ++allocations_;
    // Copy the block size in the front of returned memory.
    std::memcpy(ptr, &mem_size, sizeof(size_t));
    std::byte* offset_ptr = static_cast<std::byte*>(ptr) + SIZE_ALIGN_BYTES;
    T* data_ptr = reinterpret_cast<T*>(offset_ptr);
    return data_ptr;

  } else {

    throw std::bad_alloc{}; // required by [new.delete.single]/3

  }

}


template<class T>
void AuditMemory::deleteMem(T* ptr)
{

  if (ptr == nullptr) {

    // Freeing a null pointer is not an error.
    return;

  }

  // Recover the block size.
  size_t deallocate_size;
  std::byte* offset_ptr = reinterpret_cast<std::byte*>(ptr) - SIZE_ALIGN_BYTES;
  std::memcpy(&deallocate_size, offset_ptr, sizeof(size_t));
  ++deallocations_;
  deallocated_bytes_ += deallocate_size;
  // Thread safe in linux.
  free(static_cast<void*>(offset_ptr));

}


} // namespace.


#endif //KEL_MEM_ALLOC_H
