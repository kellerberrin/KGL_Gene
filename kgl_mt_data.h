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

#include <memory>
#include <string>
#include "kgl_logging.h"
#include "kgl_lock.h"
#include "kgl_genome_types.h"
#include "kgl_nucleotide.h"
#include "kgl_mt_numpy.h"
#include "kgl_mt_contig.h"
#include "kgl_mt_insert.h"
#include "kgl_exec_env.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Important - these typedefs define the nucleotide types being analyzed (generally DNA).
// The multi-thread locking strategy and the nucleotide data structure being updated, numpy or local.
using NumpyArray = NumpyContigMT<X86CountLock, StandardNucleotideColumn>; // Use fast asm lock and standard columns.
using LocalArray = LocalContigMT<X86CountLock, StandardNucleotideColumn>; // Use fast asm lock and standard columns.

// Define the consumer insert data structure with locking strategy.
using ConsumerInsertType = ContigInsertSequences<GranularityMutex<1000>>;

// An object to hold both the thread-safe nucleotide array and insert array for a Python Numpy.
class ConsumerNumpyRecord {

public:

  explicit ConsumerNumpyRecord(const ContigSize_t contig_size,
                               NucleotideReadCount_t *data_ptr,
                               const ContigOffset_t num_nucleotides) : nucleotide_array_(data_ptr,
                                                                                         contig_size,
                                                                                         num_nucleotides),
                                                                       insert_array_(contig_size) {}
  ~ConsumerNumpyRecord() = default;

  inline NumpyArray& getNucleotideArray() { return nucleotide_array_; }
  inline ConsumerInsertType& getInsertArray() { return insert_array_; }

private:

  NumpyArray nucleotide_array_;
  ConsumerInsertType insert_array_;

};


// Uses internal data arrays to hold the data
// An object to hold both the local thread-safe nucleotide array and thread-safe insert array.
class ConsumerLocalRecord {

public:

  explicit ConsumerLocalRecord(const ContigSize_t contig_size) : nucleotide_array_(contig_size),
                                                                 insert_array_(contig_size) {}
  ~ConsumerLocalRecord() = default;

  inline LocalArray& getNucleotideArray() { return nucleotide_array_; }
  inline ConsumerInsertType& getInsertArray() { return insert_array_; }

private:

  LocalArray nucleotide_array_;
  ConsumerInsertType insert_array_;

};

// ContigId indexed data. This is a map structure that holds the processed read data indexed by contig_id.
template<typename DataBlock> using ContigCountMap = std::map<const ContigId_t, std::shared_ptr<DataBlock>>;
template <class ContigDataBlock>
class ContigDataMap {

public:

  explicit ContigDataMap() = default;
  ~ContigDataMap() = default;

// Access functions to obtain the underlying contig block.

  ContigCountMap<ContigDataBlock>& getMap() { return contig_map_; }

  std::shared_ptr<ContigDataBlock> findContigBlock(const ContigId_t& contig_id) {

    auto result = contig_map_.find(contig_id);
    if (result != contig_map_.end()) {
      return result->second;
    }
    std::shared_ptr<ContigDataBlock> null;
    return null;

  }

  template<typename... Args>
  bool insertContig( const ContigId_t& contig_id, const ContigSize_t contig_size, Args... args) {

    std::shared_ptr<ContigDataBlock> contig_data_ptr(std::make_shared<ContigDataBlock>(contig_size, args...));
    auto result = contig_map_.insert(std::make_pair(contig_id, contig_data_ptr));
    return result.second;

  }


private:

  ContigCountMap<ContigDataBlock> contig_map_;  // Store the DNA read data for all contigs.

};

// Important - this is the data structure used in the executable.
using ContigCountData = ContigDataMap<ConsumerLocalRecord>;


}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_MT_DATA_H
