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
// Created by kellerberrin on 11/09/17.
//

#include <iostream>
#include "kgl_mt_data.h"

namespace kgl = kellerberrin::genome;

//
// Implementation of ContigMatrixMT object that provides thread safe access to the contig read count data
//

static inline void asmX86AtomicInc(kgl::NucleotideReadCount_t  *access_ptr)
{
  int inc = 1;

  __asm__ volatile("lock; xaddl %0, %1"
  : "+r" (inc), "+m" (*access_ptr) // input+output
  : // memory and condition codes changed
  : "memory", "cc"
  );

}

void kgl::ContigDataMap::addContigData( const ContigId_t& contig_id
                                      , const NucleotideReadCount_t *data_ptr
                                      , const ContigOffset_t contig_offset
                                      , const ContigOffset_t num_nucleotides)  // This is to check the numpy dimensions
{

  // Contig data matrices should be setup before threads are spawned, but let's be sure.
  std::lock_guard<std::mutex> lock(mutex_) ;

  std::unique_ptr<ContigArrayMT> contig_matrix_ptr(std::make_unique<ContigArrayMT>(log,
                                                                                   data_ptr,
                                                                                   contig_offset,
                                                                                   num_nucleotides));

  auto result = contig_map_.insert(std::make_pair(contig_id, std::move(contig_matrix_ptr)));

  if (not result.second) {

    log.error("addContigData(), Attempted to add duplicate contig; {}", contig_id);

  }

}
