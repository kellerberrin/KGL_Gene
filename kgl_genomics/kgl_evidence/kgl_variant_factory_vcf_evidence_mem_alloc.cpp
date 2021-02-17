//
// Created by kellerberrin on 15/2/21.
//

#include "kgl_variant_factory_vcf_evidence_mem_alloc.h"


namespace kgl = kellerberrin::genome;


void kgl::MemoryAllocationStrategy::allocateMemory(const MemDataUsage &mem_count, MemoryStrategy strategy) {

  // Makes the function re-entrant.
  deallocateMemory();
  strategy_ = strategy;
  // Defined Memory Strategies
  switch(strategy_) {

    case MemoryStrategy::MALLOC:
      allocateMalloc(mem_count);
      break;

    case MemoryStrategy::AUDITED_MALLOC:
      allocateAuditedMalloc(mem_count);
      break;

    case MemoryStrategy::SINGLE_MALLOC:
      allocateSingleMalloc(mem_count);
      break;

    default:
    case MemoryStrategy::AUDITED_SINGLE_MALLOC:
      allocateAuditedSingleMalloc(mem_count);
      break;

  }

}


void kgl::MemoryAllocationStrategy::deallocateMemory() {

  switch(strategy_) {

    case MemoryStrategy::MALLOC:
      deallocateMalloc();
      break;

    case MemoryStrategy::AUDITED_MALLOC:
      deallocateAuditedMalloc();
      break;

    case MemoryStrategy::SINGLE_MALLOC:
      deallocateSingleMalloc();
      break;

    default:
    case MemoryStrategy::AUDITED_SINGLE_MALLOC:
      deallocateAuditedSingleMalloc();
      break;

  }

}


// The memory allocation and deallocation strategies.
void kgl::MemoryAllocationStrategy::allocateMalloc(const MemDataUsage &mem_count) {

  // Only allocate memory if necessary.
  if (mem_count.charCount() != 0) {

    char_memory_ = new char[mem_count.charCount()];

  }
  if (mem_count.integerCount() != 0) {

    integer_memory_ = new InfoIntegerType[mem_count.integerCount()];

  }
  if (mem_count.floatCount() != 0) {

    float_memory_ = new InfoFloatType[mem_count.floatCount()];

  }
  if (mem_count.arrayCount() != 0) {

    array_memory_ = new InfoArrayIndex[mem_count.arrayCount()];

  }
  if (mem_count.stringCount() != 0) {

    string_memory_ = new std::string_view[mem_count.stringCount()];

  }

}

void kgl::MemoryAllocationStrategy::allocateAuditedMalloc(const MemDataUsage &mem_count) {

  // Only allocate memory if necessary.
  if (mem_count.charCount() != 0) {

    char_memory_ = AuditMemory::newArray<char>(mem_count.charCount());

  }
  if (mem_count.integerCount() != 0) {

    integer_memory_ = AuditMemory::newArray<InfoIntegerType>(mem_count.integerCount());

  }
  if (mem_count.floatCount() != 0) {

    float_memory_ = AuditMemory::newArray<InfoFloatType>(mem_count.floatCount()) ;

  }
  if (mem_count.arrayCount() != 0) {

    array_memory_ = AuditMemory::newArray<InfoArrayIndex>(mem_count.arrayCount());

  }
  if (mem_count.stringCount() != 0) {

    string_memory_ = AuditMemory::newArray<std::string_view>(mem_count.stringCount());

  }


}

void kgl::MemoryAllocationStrategy::allocateSingleMalloc(const MemDataUsage &mem_count) {

  const static int64_t NO_INDEX = -1;
  // Only allocate memory if necessary.
  int64_t char_index{NO_INDEX};
  int64_t integer_index{NO_INDEX};
  int64_t float_index{NO_INDEX};
  int64_t array_index{NO_INDEX};
  int64_t string_index{NO_INDEX};
  size_t mem_size{0};

  // Calculate the aligned offsets.
  if (mem_count.charCount() != 0) {

    char_index = mem_size;
    mem_size += AuditMemory::alignedArray<char>(mem_count.charCount());

  }
  if (mem_count.integerCount() != 0) {

    integer_index = mem_size;
    mem_size += AuditMemory::alignedArray<InfoIntegerType>(mem_count.integerCount());

  }
  if (mem_count.floatCount() != 0) {

    float_index = mem_size;
    mem_size += AuditMemory::alignedArray<InfoFloatType>(mem_count.floatCount());

  }
  if (mem_count.arrayCount() != 0) {

    array_index = mem_size;
    mem_size += AuditMemory::alignedArray<InfoArrayIndex>(mem_count.arrayCount());

  }
  if (mem_count.stringCount() != 0) {

    string_index = mem_size;
    mem_size += AuditMemory::alignedArray<std::string_view>(mem_count.stringCount());

  }

  // Actually allocate the memory.
  if (mem_size > 0) {

    byte_ptr_ = new std::byte[mem_size];

  }

  // Assign the data addresses using the offsets..
  if (char_index != NO_INDEX) {

    char_memory_ = reinterpret_cast<char*>(&byte_ptr_[char_index]);

  }
  if (integer_index != NO_INDEX) {

    integer_memory_ = reinterpret_cast<InfoIntegerType*>(&byte_ptr_[integer_index]);

  }
  if (float_index != NO_INDEX) {

    float_memory_ = reinterpret_cast<InfoFloatType*>(&byte_ptr_[float_index]);

  }
  if (array_index != NO_INDEX) {

    array_memory_ = reinterpret_cast<InfoArrayIndex*>(&byte_ptr_[array_index]);

  }
  if (string_index != NO_INDEX) {

    string_memory_ = reinterpret_cast<std::string_view*>(&byte_ptr_[string_index]);

  }

}

void kgl::MemoryAllocationStrategy::allocateAuditedSingleMalloc(const MemDataUsage &mem_count) {


  const static int64_t NO_INDEX = -1;
  // Only allocate memory if necessary.
  int64_t char_index{NO_INDEX};
  int64_t integer_index{NO_INDEX};
  int64_t float_index{NO_INDEX};
  int64_t array_index{NO_INDEX};
  int64_t string_index{NO_INDEX};
  size_t mem_size{0};

  // Calculate the aligned offsets.
  if (mem_count.charCount() != 0) {

    char_index = mem_size;
    mem_size += AuditMemory::alignedArray<char>(mem_count.charCount());

  }
  if (mem_count.integerCount() != 0) {

    integer_index = mem_size;
    mem_size += AuditMemory::alignedArray<InfoIntegerType>(mem_count.integerCount());

  }
  if (mem_count.floatCount() != 0) {

    float_index = mem_size;
    mem_size += AuditMemory::alignedArray<InfoFloatType>(mem_count.floatCount());

  }
  if (mem_count.arrayCount() != 0) {

    array_index = mem_size;
    mem_size += AuditMemory::alignedArray<InfoArrayIndex>(mem_count.arrayCount());

  }
  if (mem_count.stringCount() != 0) {

    string_index = mem_size;
    mem_size += AuditMemory::alignedArray<std::string_view>(mem_count.stringCount());

  }

  // Actually allocate the memory.
  if (mem_size > 0) {

    byte_ptr_ = AuditMemory::newArray<std::byte>(mem_size);

  }

  // Assign the data addresses using the offsets..
  if (char_index != NO_INDEX) {

    char_memory_ = reinterpret_cast<char*>(&byte_ptr_[char_index]);

  }
  if (integer_index != NO_INDEX) {

    integer_memory_ = reinterpret_cast<InfoIntegerType*>(&byte_ptr_[integer_index]);

  }
  if (float_index != NO_INDEX) {

    float_memory_ = reinterpret_cast<InfoFloatType*>(&byte_ptr_[float_index]);

  }
  if (array_index != NO_INDEX) {

    array_memory_ = reinterpret_cast<InfoArrayIndex*>(&byte_ptr_[array_index]);

  }
  if (string_index != NO_INDEX) {

    string_memory_ = reinterpret_cast<std::string_view*>(&byte_ptr_[string_index]);

  }

}

void kgl::MemoryAllocationStrategy::deallocateMalloc() {

  if (char_memory_ != nullptr) {

    delete[] char_memory_;
    char_memory_ = nullptr;

  }
  if (integer_memory_ != nullptr) {

    delete[] integer_memory_;
    integer_memory_ = nullptr;

  }
  if (float_memory_ != nullptr) {

    delete[] float_memory_;
    float_memory_ = nullptr;

  }
  if (array_memory_ != nullptr) {

    delete[] array_memory_;
    array_memory_ = nullptr;

  }
  if (string_memory_ != nullptr) {

    delete[] string_memory_;
    string_memory_ = nullptr;

  }

}

void kgl::MemoryAllocationStrategy::deallocateAuditedMalloc() {

  if (char_memory_ != nullptr) {

    AuditMemory::deleteArray(char_memory_);
    char_memory_ = nullptr;

  }
  if (integer_memory_ != nullptr) {

    AuditMemory::deleteArray(integer_memory_);
    integer_memory_ = nullptr;

  }
  if (float_memory_ != nullptr) {

    AuditMemory::deleteArray(float_memory_);
    float_memory_ = nullptr;

  }
  if (array_memory_ != nullptr) {

    AuditMemory::deleteArray(array_memory_);
    array_memory_ = nullptr;

  }
  if (string_memory_ != nullptr) {

    AuditMemory::deleteArray(string_memory_);
    string_memory_ = nullptr;

  }

}

void kgl::MemoryAllocationStrategy::deallocateSingleMalloc() {

  if (byte_ptr_ != nullptr) {

    delete[] byte_ptr_;
    byte_ptr_ = nullptr;

  }
  char_memory_ = nullptr;
  integer_memory_ = nullptr;
  float_memory_ = nullptr;
  array_memory_ = nullptr;
  string_memory_ = nullptr;

}

void kgl::MemoryAllocationStrategy::deallocateAuditedSingleMalloc() {

  if (byte_ptr_ != nullptr) {

    AuditMemory::deleteArray(byte_ptr_);
    byte_ptr_ = nullptr;

  }
  char_memory_ = nullptr;
  integer_memory_ = nullptr;
  float_memory_ = nullptr;
  array_memory_ = nullptr;
  string_memory_ = nullptr;

}
