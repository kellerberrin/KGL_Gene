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
#include "kgl_genome_types.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// If the target processor is X86 then use atomic incrementCountX86() with asm xaddl instruction.
// #define KGL_USE_X86_ATOMIC  // Comment out for generic mutex protected incrementCount() for all architectures.

// Only lock small sections of the data array for write access.
// This is a classic trade-off between the space taken by the mutex array and fine-grained speedup of access

class GranularityMutex {

public:

  GranularityMutex(const std::size_t array_size,const  std::size_t granularity): array_size_(array_size)
                                                                               , granularity_(granularity) {

    auto lock_count = std::size_t(array_size / granularity) + 1;
    lock_array_ = std::move(std::vector<std::mutex>{lock_count});

  }

  inline void acquire(std::size_t array_index) {

    auto idx = std::size_t(array_index / granularity_);
    lock_array_[idx].lock();

  }

  inline void release(std::size_t array_index) {

    auto idx = std::size_t(array_index / granularity_);
    lock_array_[idx].unlock();

  }

private:

  const std::size_t array_size_;
  const std::size_t granularity_;
  std::vector<std::mutex> lock_array_;

};

// This object holds the aggregated SAM read data for a DNA contiguous region. Thread safe.
// The data block is assumed to be contiguous. It is allocated from raw memory - see associated Python code.
// Important - pointer arithmetic assumes the associated numpy matrix is row-major (default).
// The data is accessed via a raw pointer - nasty but necessary. Elements are assumed set to zero.

class ContigMatrixMT {

public:

  ContigMatrixMT( Logger& logger
                , const NucleotideReadCount_t *data_ptr
                , const std::size_t rows
                , const std::size_t nucleotides  // This is to check the numpy dimensions
                , const std::size_t mutex_granularity) : log(logger)
                                                       , data_ptr_(data_ptr)
                                                       , rows_(rows)
                                                       , granularity_mutex_(rows, mutex_granularity) {

    if (nucleotides != expected_columns) {

      log.critical("Numpy array is expected to have 6 nucleotides (columns), actual columns: {}", nucleotides);
      log.critical("Serious software error - program terminates");
      std::exit(EXIT_FAILURE);

    }

    initialize(0);  // Ensure the array is properly initialized.

  }

  ~ContigMatrixMT() = default;
  ContigMatrixMT(ContigMatrixMT&&) = delete;
  ContigMatrixMT& operator=(const ContigMatrixMT&) = delete;

  void incrementCountX86(const std::size_t row, const std::size_t column);

  void incrementCount(const std::size_t row, const std::size_t column);

  inline void incrementCount(const std::size_t row, const Nucleotide_t nucleotide) {

    incrementCount(row, nucleotideToColumn(nucleotide));

  }

  const NucleotideReadCount_t  readCount(const std::size_t row, const std::size_t column) const;

  inline const NucleotideReadCount_t  readCount(const std::size_t row, const Nucleotide_t nucleotide) const {

    return readCount(row, nucleotideToColumn(nucleotide));

  }

  inline const ContigSize_t contigSize() const { return rows_; }

private:

  void initialize(const NucleotideReadCount_t initial_value);

  std::size_t nucleotideToColumn(const Nucleotide_t nucleotide) const;

  static constexpr std::size_t expected_columns = 6; // A, C, G, T/U, -, + (in that order). See nucleotideToColumn(const char).

  // Yes folks, that's a raw pointer (watch it snort and shake).
  // Never delete[] or use to initialize a smart pointer (same thing).
  // The data actually belongs to a numpy passed in from Python.
  const NucleotideReadCount_t  *data_ptr_; // Snort, shake, snort, shake ...
  const std::size_t rows_;      // Number of of nucleotides in the contig.
  mutable std::mutex mutex_;  // Used to initialize the data.
  GranularityMutex granularity_mutex_; // Fine grained thread access to the underlying data structure.

  Logger& log;

};

class ContigDataMap {

public:

  explicit ContigDataMap(Logger& logger) : log(logger) {}
  ~ContigDataMap() = default;
  ContigDataMap(ContigMatrixMT&&) = delete;
  ContigDataMap& operator=(const ContigMatrixMT&) = delete;

  void contigIncrementCount( const ContigId_t& contig_id
                           , const ContigOffset_t contig_offset
                           , const Nucleotide_t nucleotide);

  void addContigData( const ContigId_t& contig_id
                    , const NucleotideReadCount_t *data_ptr
                    , const ContigOffset_t contig_offset
                    , const ContigOffset_t num_nucleotides);  // This is to check the numpy dimensions

  ContigMatrixMT& getContig(const ContigId_t& contig_id); // Access function.

private:

  static constexpr std::size_t lock_granularity = 1000;

  std::map<const ContigId_t, std::unique_ptr<ContigMatrixMT>> contig_map_;  // Store the DNA read data for all contigs.
  mutable std::mutex mutex_;  // Used to add contigs.

  Logger& log;

};


// Class to hold enqueued nucleotide sequences to be inserted in the genome model.
// The consumer threads enqueue inserted sequences along with contig id and offset.
// When the consumer threads complete, the mainline process transfers this information back to the Python
// program so that it can construct a data array of insertions for each contig.

class InsertQueue {

public:

  InsertQueue() = default;
  ~InsertQueue() = default;
  InsertQueue(const InsertQueue&) = delete;
  InsertQueue(InsertQueue&&) = delete;
  InsertQueue& operator=(const InsertQueue&) = delete;

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
