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
// Created by kellerberrin on 10/09/17.
//

#ifndef KGL_MT_DATA_H
#define KGL_MT_DATA_H

#include <cstdint>
#include <memory>
#include <string>
#include <queue>
#include "kgl_logging.h"
#include "kgl_lock.h"
#include "kgl_genome_types.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Standard contig column layout. Implement the nucleotides as "A", "C", "G", "T"/"U", "-", "+" in that order.
class StandardNucleotideColumn {

public:

  StandardNucleotideColumn(Logger& logger) : log(logger) {}

  static constexpr ContigOffset_t nucleotides = 7;
  static constexpr Nucleotide_t delete_nucleotide = '-';
  static constexpr Nucleotide_t insert_sequence = '+';

  ContigOffset_t nucleotideToColumn(const Nucleotide_t nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case 'A':
      case 'a': return 0;

      case 'C':
      case 'c': return 1;

      case 'G':
      case 'g': return 2;

      case 'U':
      case 'u':
      case 'T':
      case 't': return 3;

      case 'N':
      case 'n': return 4;

      case '-': return 5;

      case '+': return 6;

      default:
        log.critical("nucleotideToColumn(), Count data array accessed with unknown nucleotide: {}", nucleotide);
        return 0; // Never reached, to keep the compiler happy.

    }

  }

private:

  Logger& log;

};

// This object holds the aggregated SAM read data for a DNA contiguous region. Thread safe LockStrategy.
// The data block is assumed to be contiguous. It is allocated from raw memory - see associated Python code.
// Important - pointer arithmetic assumes the associated numpy matrix is row-major (default).
// The data is accessed via a raw pointer - nasty but necessary. Elements are assumed set to zero.

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

    if (nucleotides != NucleotideColumn::nucleotides) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::nucleotides;
      log.critical("Numpy array is expected to have: {} nucleotides (columns), columns specified: {}"
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
    if (column >= NucleotideColumn::nucleotides) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::nucleotides;
      log.critical("Invalid access in incrementCount(); Nucleotide column index: {} >= Number Nucleotides: {}"
                   , column, nucleotide_columns);

    }

    // pointer arithmetic, stride is expected_columns
    NucleotideReadCount_t* access_ptr = const_cast<NucleotideReadCount_t *>(data_ptr_)
                                        + ((contig_offset * NucleotideColumn::nucleotides) + column);

    lock_strategy_.incrementCount(*access_ptr, contig_offset);

  }

  // Reads are not locked. Do not call while updating read counts.
  inline const NucleotideReadCount_t readCount( const ContigOffset_t contig_offset
                                              , const Nucleotide_t nucleotide) const {

    std::size_t column = nucleotide_column_.nucleotideToColumn(nucleotide);

    if (contig_offset >= contig_size_) {

      log.critical("Invalid access in incrementCount(); Contig index: {} >= Contig size: {}"
                  , contig_offset, contig_size_);

    }
    if (column >= NucleotideColumn::nucleotides) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::nucleotides;
      log.critical("Invalid access in incrementCount(); Nucleotide column index: {} >= Number Nucleotides: {}"
                  , column, nucleotide_columns);

    }

    // pointer arithmetic, stride is number of nucleotides
    const NucleotideReadCount_t  *access_ptr = data_ptr_ + ((contig_offset * NucleotideColumn::nucleotides) + column);
    return *access_ptr;  // read access - no mutex.

  }

  inline const ContigSize_t contigSize() const { return contig_size_; }

private:

  Logger& log;

  NucleotideColumn nucleotide_column_;
  // Yes folks, that's a raw pointer (watch it snort and shake).
  // The data actually belongs to a numpy passed in from Python.
  const NucleotideReadCount_t *data_ptr_; // Snort, shake, snort, shake ...
  const std::size_t contig_size_;      // Number of of nucleotides in the contig.
  mutable std::mutex mutex_;  // Used to initialize the data.
  LockStrategy lock_strategy_; // Fine grained thread access to the underlying data structure (or X86 asm).

  void initialize(const NucleotideReadCount_t initial_value) {

    // Initialization happens before threads are spawned but we lock the data just to make sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    for (ContigOffset_t contig_offset = 0; contig_offset < contig_size_; ++contig_offset) {
      for (ContigOffset_t column = 0; column < NucleotideColumn::nucleotides; ++column) {
        // pointer arithmetic, stride is expected_columns
        NucleotideReadCount_t  *access_ptr = const_cast<NucleotideReadCount_t  *>(data_ptr_)
                                             + ((contig_offset * NucleotideColumn::nucleotides) + column);
        (*access_ptr) = initial_value;
      }

    }

  }

};

// Contig Id indexed Contig Arrays.

using ContigArrayMT = NumpyContigMT<NullMutex<NucleotideReadCount_t >, StandardNucleotideColumn>; // Use fast asm Mutex and standard columns.
using ContigMap = std::map<const ContigId_t, std::unique_ptr<ContigArrayMT>>;

class ContigDataMap {

public:

  explicit ContigDataMap(Logger& logger) : log(logger) {}
  ~ContigDataMap() = default;

  void addContigData( const ContigId_t& contig_id
                    , const NucleotideReadCount_t *data_ptr
                    , const ContigOffset_t contig_offset
                    , const ContigOffset_t num_nucleotides);  // This is to check the numpy dimensions

// Access function to obtain the underlying contig block.
  inline ContigMap::iterator getContig( const ContigId_t& contig_id) { return contig_map_.find(contig_id); }
  inline ContigMap::const_iterator notFound() { return contig_map_.end(); }
  inline ContigArrayMT& getMatrix(ContigMap::iterator& map_ptr) { return *(map_ptr->second); }

private:

  Logger& log;

  ContigMap contig_map_;  // Store the DNA read data for all contigs.
  mutable std::mutex mutex_;  // Used to add contigs.

};

// Class to hold enqueued nucleotide sequences to be inserted in the genome model.
// The consumer threads enqueue inserted sequences along with contig id and offset.
// When the consumer threads complete, the mainline process transfers this information back to the Python
// program so that it can construct a data array of insertions for each contig.

class InsertQueue {

public:

  InsertQueue() = default;
  ~InsertQueue() = default;

  // Thread safe
  void push(const ContigId_t& contig_id, const ContigOffset_t& contig_offset, const Sequence_t& sequence) {

    std::lock_guard<std::mutex> lock(mutex_);

    contig_id_vector_.push_back(contig_id);
    contig_offset_vector_.push_back(contig_offset);
    sequence_vector_.push_back(sequence);

  }

  // Thread unsafe. Access to this data structure may only be performed after the spawned threads have joined.
  const std::vector<std::string>& getQueueContigs() {

    return contig_id_vector_;

  }

  const std::vector<std::size_t>& getQueueOffsets() {

    return contig_offset_vector_;

  }

  const std::vector<std::string>& getQueueSequences() {

    return sequence_vector_;

  }

private:

  mutable std::mutex mutex_;

  std::vector<ContigId_t> contig_id_vector_;
  std::vector<ContigOffset_t> contig_offset_vector_;
  std::vector<Sequence_t> sequence_vector_;

};

}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_MT_DATA_H
