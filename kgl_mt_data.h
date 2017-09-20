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
#include "kgl_nucleotide.h"
#include "kgl_mt_numpy.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Important - these typedefs define the nucleotide types being analyzed (generally DNA).
// The multi-thread locking strategy and the data structure being updated, numpy or local.

using ContigArrayMT = NumpyContigMT<X86Mutex, StandardNucleotideColumn>; // Use fast asm Mutex and standard columns.
using ContigMap = std::map<const ContigId_t, std::unique_ptr<ContigArrayMT>>;


// ContigId indexed ContigArrayMT.
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
