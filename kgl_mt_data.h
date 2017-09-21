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
#include <iostream>
#include <memory>
#include <string>
#include <queue>
#include "kgl_logging.h"
#include "kgl_lock.h"
#include "kgl_genome_types.h"
#include "kgl_nucleotide.h"
#include "kgl_mt_numpy.h"
#include "kgl_mt_contig.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// The insert sequence data structure is an contig sized array of maps indexed by contig offset.
using OffsetSequenceMap = std::map<const Sequence_t, NucleotideReadCount_t>;

template <class LockingStrategy = ArrayMutex>    // Specify a GranularityMutex or an ArrayMutex.
class ContigInsertSequences {

public:

  ContigInsertSequences(Logger& logger, ContigSize_t contig_size)
      : log(logger),
        insert_array_ptr_(std::make_unique<OffsetSequenceMap[]>(contig_size)),
                                                    mutex_array_(contig_size),
                                                    contig_size_(contig_size) {}
  ~ContigInsertSequences() = default;

  inline void insertSequence(ContigOffset_t contig_offset, const Sequence_t& sequence) {


    if (contig_offset >= contig_size_) {

      log.error("Insert offset: {} exceeds contig size: {}", contig_offset, contig_size_);

    }

    // Lock the data structure.
    mutex_array_.acquire(contig_offset);

    // Check if the inserted sequence already exists.
    auto found = insert_array_ptr_[contig_offset].find(sequence);

    // Add the sequence if not found.
    if (found == insert_array_ptr_[contig_offset].end()) {

      auto result = insert_array_ptr_[contig_offset].insert(std::make_pair(sequence, 1));

      if (not result.second) {

        log.error("insertSequence(), Attempted to add duplicate sequence: {}, at offset: {}", sequence, contig_offset);
      }

    } else { // Sequence was found, so increment counter.

      ++(found->second);

    }

    mutex_array_.release(contig_offset);

  }


private:

  Logger& log;

  std::unique_ptr<OffsetSequenceMap[]> insert_array_ptr_;   // points to contig sized array of sequence maps.
  mutable LockingStrategy mutex_array_;           // Locks the data structure for update.
  ContigSize_t contig_size_;


};

// ContigId indexed data.
template <class ContigDataBlock> class ContigDataMap {

public:

  using ContigMap = std::map<const ContigId_t, std::unique_ptr<ContigDataBlock>>;

  explicit ContigDataMap(Logger& logger) : log(logger) {}
  ~ContigDataMap() = default;


  void addContigBlock( const ContigId_t& contig_id, std::unique_ptr<ContigDataBlock>& data_ptr)
  {

    // Contig data matrices should be setup before threads are spawned, but let's be sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    auto result = contig_map_.insert(std::make_pair(contig_id, std::move(data_ptr)));

    if (not result.second) {

      log.error("addContigData(), Attempted to add duplicate contig; {}", contig_id);

    }

  }

// Access function to obtain the underlying contig block.
  inline typename ContigMap::iterator getContig( const ContigId_t& contig_id) { return contig_map_.find(contig_id); }
  inline typename ContigMap::const_iterator notFound() { return contig_map_.end(); }
  inline ContigDataBlock& getMatrix(typename ContigMap::iterator& map_ptr) { return *(map_ptr->second); }

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
