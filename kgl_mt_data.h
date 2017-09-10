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
#include "kgl_logging.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

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
                , const uint32_t *data_ptr
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

  inline void incrementCount(const std::size_t row, const std::size_t column) {

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
    uint32_t *access_ptr = const_cast<uint32_t *>(data_ptr_) + ((row * expected_columns) + column);
    granularity_mutex_.acquire(row);  // Enforce thread protection.
    ++(*access_ptr);
    granularity_mutex_.release(row);

  }

  inline void incrementCount(const std::size_t row, const char nucleotide) {

    incrementCount(row, nucleotideToColumn(nucleotide));

  }


  inline const uint32_t readCount(const std::size_t row, const std::size_t column) const {

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
    const uint32_t *access_ptr = data_ptr_ + ((row * expected_columns) + column);
    return *access_ptr;  // read access - no mutex.

  }

  inline const uint32_t readCount(const std::size_t row, const char nucleotide) const {

    return readCount(row, nucleotideToColumn(nucleotide));

  }

private:

  void initialize(const uint32_t initial_value) {

    // Initialization happens before threads are spawned but we lock the data just to make sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    for (std::size_t row = 0; row < rows_; ++row) {
      for (std::size_t column = 0; column < expected_columns; ++column) {
        // pointer arithmetic, stride is expected_columns
        uint32_t *access_ptr = const_cast<uint32_t *>(data_ptr_) + ((row * expected_columns) + column);
        (*access_ptr) = initial_value;
      }

    }

  }

  inline std::size_t nucleotideToColumn(const char nucleotide) const {

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

      case '-':
        column = 4;
        break;

      case '+':
        column = 5;
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

  static constexpr std::size_t expected_columns = 6; // A, C, G, T/U, -, + (in that order). See nucleotideToColumn(const char).

  // Yes folks, that's a raw pointer (watch it snort and shake).
  // Never delete[] or use to initialize a smart pointer (same thing).
  // The data actually belongs to a numpy passed in from Python.
  const uint32_t *data_ptr_;
  const std::size_t rows_;
  mutable std::mutex mutex_;  // Used to initialize the data.
  GranularityMutex granularity_mutex_;

  Logger& log;

};

class ContigDataMap {

public:

  explicit ContigDataMap(Logger& logger) : log(logger) {}
  ~ContigDataMap() = default;
  ContigDataMap(ContigMatrixMT&&) = delete;
  ContigDataMap& operator=(const ContigMatrixMT&) = delete;

  void contigIncrementCount(const std::string& contig_id, const std::size_t row, const char nucleotide) {

    auto search = contig_map_.find(contig_id);

    if(search != contig_map_.end()) {

      search->second.incrementCount(row, nucleotide);

    }
    else {

      log.error("contigIncrementCount(), No data array exists for contig; {}", contig_id);
      log.error("Program exists");
      std::exit(EXIT_FAILURE);

    }

  }

  void addContigData( const std::string& contig_id
                    , const uint32_t *data_ptr
                    , const std::size_t rows
                    , const std::size_t nucleotides)  // This is to check the numpy dimensions
  {

    // Contig data matrices should be setup before threads are spawned, but let's be sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    using key_val = std::pair<const std::string, ContigMatrixMT>;
    ContigMatrixMT contig_matrix(log, data_ptr, rows, nucleotides, granularity);
    key_val kv_pair = ( contig_id, contig_matrix);
    auto result = contig_map_.insert(kv_pair);

    if (not result.second) {

      log.error("addContigData(), Attempted to add duplicate contig; {}", contig_id);
      log.error("Program exists");
      std::exit(EXIT_FAILURE);

    }
  }


private:

  static constexpr std::size_t granularity = 1000;

  std::map<const std::string, ContigMatrixMT> contig_map_;  // Store the DNA read data for all contigs.
  mutable std::mutex mutex_;  // Used to add contigs.

  Logger& log;

};


}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_MT_DATA_H
