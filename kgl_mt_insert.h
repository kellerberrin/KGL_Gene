//
// Created by kellerberrin on 11/10/17.
//

#ifndef KGL_MT_INSERT_H
#define KGL_MT_INSERT_H

#include <memory>
#include <string>
#include <vector>
#include "kgl_lock.h"
#include "kgl_genome_types.h"
#include "kgl_alphabet_extend.h"
#include "kgl_exec_env.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Class to hold inserted nucleotide sequences to be passed back to Python.
// This object is not updated directly from processed SAM records, but is populated
// from the ContigInsertSequences class below.
class InsertQueue {

public:

  InsertQueue() = default;
  ~InsertQueue() = default;

  // Thread safe
  void push(const ContigId_t& contig_id, const ContigOffset_t& contig_offset, const AlphabetSequence_t& sequence) {

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
  std::vector<AlphabetSequence_t> sequence_vector_;

};


// The insert sequence data structure is an contig sized array of maps indexed by contig offset.
using OffsetSequenceMap = std::map<const AlphabetSequence_t, NucleotideReadCount_t>;

template <class LockingStrategy = ArrayMutex>    // Specify a GranularityMutex or an ArrayMutex.
class ContigInsertSequences {

public:

  ContigInsertSequences(ContigSize_t contig_size)
      :  insert_array_ptr_(std::make_unique<OffsetSequenceMap[]>(contig_size)),
         mutex_array_(contig_size),
         contig_size_(contig_size) {}
  ~ContigInsertSequences() = default;

  inline void insertSequence(ContigOffset_t contig_offset, const AlphabetSequence_t& sequence) {

    if (contig_offset >= contig_size_) {

      ExecEnv::log().error("Insert offset: {} exceeds contig size: {}", contig_offset, contig_size_);

    }

    // Lock the data structure.
    mutex_array_.acquire(contig_offset);

    // Check if the inserted sequence already exists.
    auto found = insert_array_ptr_[contig_offset].find(sequence);

    // Add the sequence if not found.
    if (found == insert_array_ptr_[contig_offset].end()) {

      auto result = insert_array_ptr_[contig_offset].insert(std::make_pair(sequence, 1));

      if (not result.second) {

        ExecEnv::log().error("insertSequence(), Attempted to add duplicate sequence: {}, at offset: {}",
                             sequence, contig_offset);
      }

    } else { // Sequence was found, so increment counter.

      ++(found->second);

    }

    mutex_array_.release(contig_offset);

  }

  inline void convertToQueue(const ContigId_t& contig_id, InsertQueue& insert_queue) {

    for (ContigOffset_t offset = 0; offset < contig_size_; ++offset) {

      for (auto sequence : insert_array_ptr_[offset]) {

        for (NucleotideReadCount_t count = 0; count < sequence.second; count++) {

          insert_queue.push(contig_id, offset, sequence.first);

        }

      }

    }

  }

private:

  std::unique_ptr<OffsetSequenceMap[]> insert_array_ptr_;   // points to contig sized array of sequence maps.
  mutable LockingStrategy mutex_array_;           // Locks the data structure for update.
  ContigSize_t contig_size_;

};


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_MT_INSERT_H
