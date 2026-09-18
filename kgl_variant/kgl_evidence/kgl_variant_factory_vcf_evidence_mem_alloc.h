//
// Created by kellerberrin on 15/2/21.
//

#ifndef KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEM_ALLOC_H
#define KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEM_ALLOC_H

#include "kel_mem_alloc.h"
#include "kgl_variant_factory_vcf_evidence_data.h"
#include "kgl_variant_factory_vcf_evidence_memory.h"



namespace kellerberrin::genome {   //  organization level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implements the memory allocation strategy of the evidence object.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Defined Memory Strategies
enum class MemoryStrategy {  MALLOC, AUDITED_MALLOC, SINGLE_MALLOC, AUDITED_SINGLE_MALLOC };


class MemoryAllocationStrategy {

public:

  MemoryAllocationStrategy() = default;
  ~MemoryAllocationStrategy() { deallocateMemory(); }

  [[nodiscard]] char* charMemory() const { return char_memory_; }
  [[nodiscard]] InfoIntegerType* integerMemory() const { return integer_memory_; }
  [[nodiscard]] InfoFloatType* floatMemory() const { return float_memory_; };
  [[nodiscard]] InfoArrayIndex* arrayMemory() const { return array_memory_; }
  [[nodiscard]] std::string_view* stringMemory() const { return string_memory_; }
  [[nodiscard]] MemoryStrategy memoryStrategy() { return strategy_; }

  void allocateMemory(const MemDataUsage &mem_count, MemoryStrategy strategy);
  void deallocateMemory();

private:

  char* char_memory_{nullptr};
  InfoIntegerType* integer_memory_{nullptr};
  InfoFloatType* float_memory_{nullptr};
  InfoArrayIndex* array_memory_{nullptr};
  std::string_view* string_memory_{nullptr};
  std::byte* byte_ptr_{nullptr};
  MemoryStrategy strategy_{MemoryStrategy::AUDITED_MALLOC};

  // The allocation and deallocation strategies.
  void allocateMalloc(const MemDataUsage &mem_count);
  void allocateAuditedMalloc(const MemDataUsage &mem_count);
  void allocateSingleMalloc(const MemDataUsage &mem_count);
  void allocateAuditedSingleMalloc(const MemDataUsage &mem_count);

  void deallocateMalloc();
  void deallocateAuditedMalloc();
  void deallocateSingleMalloc();
  void deallocateAuditedSingleMalloc();


};



} //namespace

#endif //KGL_VARIANT_FACTORY_VCF_EVIDENCE_MEM_ALLOC_H
