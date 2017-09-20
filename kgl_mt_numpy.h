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
// Created by kellerberrin on 20/09/17.
//

#ifndef KGL_MT_NUMPY_H
#define KGL_MT_NUMPY_H

#include <cstdint>
#include <memory>
#include <string>
#include <queue>
#include "kgl_logging.h"
#include "kgl_lock.h"
#include "kgl_genome_types.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// This object holds the aggregated SAM read data for a DNA contiguous memory region passed in by a Python
// numpy.  Important - The numpy array is assumed to be a contiguous block of unsigned long (NucleotideReadCount_t).
// The equivalent numpy type is 'np.dtype("uint32")' (see the associated Python code).
// The pointer arithmetic assumes the associated numpy matrix is row-major (numpy default) and is contig_size rows by
// NucleotideColumns:: NUCLEOTIDE_COLUMNS columns passed in as a template argument. If any of these these conditions are
// not observed, then a segmentation fault and sorrow will surely follow.
// The data is accessed via a raw pointer - nasty but necessary. The multi-thread locking strategy is passed in as
// a template.

template <class LockStrategy, class NucleotideColumn>
class NumpyContigMT {

public:

  NumpyContigMT(  Logger& logger
      , const NucleotideReadCount_t *data_ptr
      , const ContigSize_t contig_size
      , const std::size_t nucleotides) :  log(logger),
                                          nucleotide_column_(log),
                                          data_ptr_(data_ptr),
                                          contig_size_(contig_size),
                                          lock_strategy_(contig_size) {

    if (nucleotides != NucleotideColumn::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::NUCLEOTIDE_COLUMNS;
      log.critical("Numpy array is expected to have: {} NUCLEOTIDE_COLUMNS (columns), columns specified: {}"
          , nucleotide_columns, nucleotides);

    }

    initialize(0);  // Ensure the array is properly initialized.

  }
  ~NumpyContigMT() = default;

  // Threadsafe read count increment.
  inline void incrementCount(const ContigOffset_t contig_offset, const Nucleotide_t nucleotide) {

    std::size_t column =  nucleotide_column_.nucleotideToColumn(nucleotide);

    if (contig_offset >= contig_size_) {

      log.critical("Invalid access in incrementCount(); Contig index: {} >= Contig size: {}"
          , contig_offset, contig_size_);

    }
    if (column >= NucleotideColumn::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::NUCLEOTIDE_COLUMNS;
      log.critical("Invalid access in incrementCount(); Nucleotide column index: {} >= Number Nucleotides: {}"
          , column, nucleotide_columns);

    }

    // pointer arithmetic, stride is expected_columns
    NucleotideReadCount_t* access_ptr = const_cast<NucleotideReadCount_t *>(data_ptr_)
                                        + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);

    lock_strategy_.incrementCount(*access_ptr, contig_offset);

  }

  inline void incrementDelete(const ContigOffset_t contig_offset) {

    incrementCount(contig_offset, NucleotideColumn::DELETE_NUCLEOTIDE);

  }

  inline void incrementInsert(const ContigOffset_t contig_offset) {

    incrementCount(contig_offset, NucleotideColumn::INSERT_SEQUENCE );

  }

  // Reads are not locked. Do not call while updating read counts.
  inline const NucleotideReadCount_t readCount( const ContigOffset_t contig_offset,
                                                const Nucleotide_t nucleotide) const {

    std::size_t column = nucleotide_column_.nucleotideToColumn(nucleotide);

    if (contig_offset >= contig_size_) {

      log.critical("Invalid access in incrementCount(); Contig index: {} >= Contig size: {}"
          , contig_offset, contig_size_);

    }
    if (column >= NucleotideColumn::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::NUCLEOTIDE_COLUMNS;
      log.critical("Invalid access in incrementCount(); Nucleotide column index: {} >= Number Nucleotides: {}"
          , column, nucleotide_columns);

    }

    // pointer arithmetic, stride is number of NUCLEOTIDE_COLUMNS
    const NucleotideReadCount_t  *access_ptr = data_ptr_ + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);
    return *access_ptr;  // read access - no mutex.

  }

  inline const ContigSize_t contigSize() const { return contig_size_; }

private:

  Logger& log;

  NucleotideColumn nucleotide_column_;
  // Yes folks, that's a raw pointer (watch it snort and shake).
  // The data actually belongs to a numpy passed in from Python.
  const NucleotideReadCount_t *data_ptr_; // Snort, shake, snort, shake ...
  const ContigSize_t contig_size_;      // Number of of NUCLEOTIDE_COLUMNS in the contig.
  mutable std::mutex mutex_;  // Used to initialize the data.
  LockStrategy lock_strategy_; // The template multi-thread locking strategy.

  void initialize(const NucleotideReadCount_t initial_value) {

    // Initialization happens before threads are spawned but we lock the data just to make sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    for (ContigOffset_t contig_offset = 0; contig_offset < contig_size_; ++contig_offset) {
      for (ContigOffset_t column = 0; column < NucleotideColumn::NUCLEOTIDE_COLUMNS; ++column) {
        // pointer arithmetic, stride is expected_columns
        NucleotideReadCount_t  *access_ptr = const_cast<NucleotideReadCount_t  *>(data_ptr_)
                                             + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);
        (*access_ptr) = initial_value;
      }

    }

  }

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_MT_NUMPY_H
