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
#include "kgl_alphabet_extend.h"
#include "kgl_mt_numpy.h"
#include "kgl_mt_contig.h"
#include "kgl_mt_insert.h"
#include "kgl_exec_env.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Important - these typedefs define the nucleotide types being analyzed (generally DNA).
// The multi-thread locking strategy and the nucleotide data structure being updated, numpy or local.
using NumpyArray = NumpyContigMT<X86CountLock, ExtendDNA5>; // Use fast asm lock and standard columns.
using LocalArray = LocalContigMT<X86CountLock, ExtendDNA5>; // Use fast asm lock and standard columns.

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

  inline const NumpyArray& getNucleotideArray() const { return nucleotide_array_; }
  inline NumpyArray& getNucleotideArray() { return nucleotide_array_; }
  inline const ConsumerInsertType& getInsertArray() const { return insert_array_; }
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

  inline const LocalArray& getNucleotideArray() const { return nucleotide_array_; }
  inline LocalArray& getNucleotideArray() { return nucleotide_array_; }
  inline const ConsumerInsertType& getInsertArray() const { return insert_array_; }
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

// File name
  void fileName(const std::string& file_name) { file_name_ = file_name;}
  const std::string& fileName() const { return file_name_; }

// Access functions to obtain the underlying contig block.

  const ContigCountMap<ContigDataBlock>& getMap() const { return contig_map_; }

  std::shared_ptr<ContigDataBlock> findContigBlock(const ContigId_t& contig_id) {

    auto result = contig_map_.find(contig_id);
    if (result != contig_map_.end()) {
      return result->second;
    }
    std::shared_ptr<ContigDataBlock> null;
    return null;

  }

  std::shared_ptr<const ContigDataBlock> findContigBlock(const ContigId_t& contig_id) const {

    auto result = contig_map_.find(contig_id);
    if (result != contig_map_.end()) {
      return result->second;
    }
    std::shared_ptr<const ContigDataBlock> null;
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
  std::string file_name_;  // File used to create the read data

};

// Important - this is the data structure used in the executable.
using ContigCountData = ContigDataMap<ConsumerLocalRecord>;


}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_MT_DATA_H
