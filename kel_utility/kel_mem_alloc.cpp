//
// Created by kellerberrin on 15/2/21.
//

#include "kel_exec_env.h"
#include "kel_mem_alloc.h"

namespace kel = kellerberrin;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Align the size of the requested memory block with alignof(max_align_t) by increasing
// the size of the requested block so that the end of the block has memory alignment alignof(max_align_t).
// This means if we concatenate requested aligned memory blocks the
// start of the next block will also be on a alignof(max_align_t) boundary.
// Thus we can take the address of the concatenated block as if it were
// a separate malloc (or ::new) since these functions always return a pointer
// with memory alignment alignof(max_align_t).
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

size_t kel::AuditMemory::alignedSize(size_t mem_size) {

  size_t alignment = alignof(std::max_align_t);
  return ((mem_size + alignment - 1) / alignment) * alignment;

}

