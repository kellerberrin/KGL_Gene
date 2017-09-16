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
// Implementation of ContigMatrixMT object that provides thread safe access to the contig Python numpy.
//

static inline void asmX86AtomicInc(kgl::NucleotideReadCount_t  *access_ptr)
{
  int inc = 1;

  __asm__ volatile("lock; xaddl %0, %1"
  : "=r" (inc), "+m" (*access_ptr) // input+output
  : // memory and condition codes changed
  : "memory", "cc"
  );

}

void kgl::ContigMatrixMT::incrementCountX86(const std::size_t row, const std::size_t column) {


  if (row >= rows_) {
    log.error("Invalid access in incrementCount(); Row index: {} >= Array size: {}", row, rows_);
    log.error("Program exits");
    std::exit(EXIT_FAILURE);
  }
  if (column >= expected_columns) {
    log.error("Invalid access in incrementCount(); Column index: {} >= Number Columns: {}", column, expected_columns);
    log.error("Program exits");
    std::exit(EXIT_FAILURE);
  }

  // pointer arithmetic, stride is expected_columns
  NucleotideReadCount_t  *access_ptr = const_cast<NucleotideReadCount_t *>(data_ptr_) + ((row * expected_columns) + column);
  asmX86AtomicInc(access_ptr);

}

void kgl::ContigMatrixMT::incrementCount(const std::size_t row, const std::size_t column) {

  if (row >= rows_) {
    log.error("Invalid access in incrementCount(); Row index: {} >= Array size: {}", row, rows_);
    log.error("Program exits");
    std::exit(EXIT_FAILURE);
  }
  if (column >= expected_columns) {
    log.error("Invalid access in incrementCount(); Column index: {} >= Number Columns: {}", column, expected_columns);
    log.error("Program exits");
    std::exit(EXIT_FAILURE);
  }

  // pointer arithmetic, stride is expected_columns
  NucleotideReadCount_t  *access_ptr = const_cast<NucleotideReadCount_t *>(data_ptr_) + ((row * expected_columns) + column);
  granularity_mutex_.acquire(row);  // Enforce thread protection.
  ++(*access_ptr);
  granularity_mutex_.release(row);

}

const kgl::NucleotideReadCount_t  kgl::ContigMatrixMT::readCount(const std::size_t row, const std::size_t column) const {

  if (row >= rows_) {
    log.error("Invalid access in readCount(); Row index: {} >= Array size: {}", row, rows_);
    log.error("Program exists");
    std::exit(EXIT_FAILURE);
  }
  if (column >= expected_columns) {
    log.error("Invalid access in readCount(); Column index: {} >= Number Columns: {}", column, expected_columns);
    log.error("Program exists");
    std::exit(EXIT_FAILURE);
  }

  // pointer arithmetic, stride is expected_columns
  const NucleotideReadCount_t  *access_ptr = data_ptr_ + ((row * expected_columns) + column);
  return *access_ptr;  // read access - no mutex.

}

void kgl::ContigMatrixMT::initialize(const NucleotideReadCount_t initial_value) {

  // Initialization happens before threads are spawned but we lock the data just to make sure.
  std::lock_guard<std::mutex> lock(mutex_) ;

  for (std::size_t row = 0; row < rows_; ++row) {
    for (std::size_t column = 0; column < expected_columns; ++column) {
      // pointer arithmetic, stride is expected_columns
      NucleotideReadCount_t  *access_ptr = const_cast<NucleotideReadCount_t  *>(data_ptr_) + ((row * expected_columns) + column);
      (*access_ptr) = initial_value;
    }

  }

}

std::size_t kgl::ContigMatrixMT::nucleotideToColumn(const Nucleotide_t nucleotide) const {

  std::size_t column;

  // Translate the nucleotide to an array column
  switch (nucleotide) {

    case 'A':
    case 'a':
      column = 0;
      break;

    case 'C':
    case 'c':
      column = 1;
      break;

    case 'G':
    case 'g':
      column = 2;
      break;

    case 'T':
    case 't':
      column = 3;
      break;

    case 'N':
    case 'n':
      column = 4;
      break;

    case '-':
      column = 5;
      break;

    case '+':
      column = 6;
      break;

    case 'U':
    case 'u':
      column = 3;
      break;

    default:
      log.error("nucleotideToColumn(), Count data array accessed with unknown nucleotide: {}", nucleotide);
      log.error("Does the Fasta file match the SAM file alignment?");
      log.error("Program exists");
      std::exit(EXIT_FAILURE);

  }

  return column;

}

//
// Implements the ContigDataMap class which provides contig indexed data arrays (underlying numpys).
//

constexpr std::size_t kgl::ContigDataMap::lock_granularity;  // Only lock small sections of the read data matrix.

// Access function to obtain the underlying contig block.
kgl::ContigMatrixMT& kgl::ContigDataMap::getContig( const ContigId_t& contig_id) {

  auto search = contig_map_.find(contig_id);

  if(search != contig_map_.end()) {

    return *(search->second);

  }
  else {

    log.error("contigIncrementCount(), No data array exists for contig; {}", contig_id);
    log.error("Program exists");
    std::exit(EXIT_FAILURE);

  }

}

void kgl::ContigDataMap::contigIncrementCount( const ContigId_t& contig_id
                                             , const ContigOffset_t contig_offset
                                             , const Nucleotide_t nucleotide) {


#ifdef KGL_USE_X86_ATOMIC
  getContig(contig_id).incrementCountX86(contig_offset, nucleotide);
#else
  getContig(contig_id).incrementCount(contig_offset, nucleotide);
#endif

}

void kgl::ContigDataMap::addContigData( const ContigId_t& contig_id
                                      , const NucleotideReadCount_t *data_ptr
                                      , const ContigOffset_t contig_offset
                                      , const ContigOffset_t num_nucleotides)  // This is to check the numpy dimensions
{

  // Contig data matrices should be setup before threads are spawned, but let's be sure.
  std::lock_guard<std::mutex> lock(mutex_) ;

  std::unique_ptr<ContigMatrixMT> contig_matrix_ptr(std::make_unique<ContigMatrixMT>(log,
                                                                                     data_ptr,
                                                                                     contig_offset,
                                                                                     num_nucleotides,
                                                                                     lock_granularity));

  auto result = contig_map_.insert(std::make_pair(contig_id, std::move(contig_matrix_ptr)));

  if (not result.second) {

    log.error("addContigData(), Attempted to add duplicate contig; {}", contig_id);
    log.error("Program exists");
    std::exit(EXIT_FAILURE);

  }

}
