//
// Created by kellerberrin on 20/09/17.
//

#ifndef KGL_MT_CONTIG_H
#define KGL_MT_CONTIG_H

#include <iostream>
#include <cstdint>
#include <memory>
#include <string>
#include <queue>
#include "kgl_logging.h"
#include "kgl_mt_numpy.h"
#include "kgl_genome_types.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


template <class LockStrategy, class NucleotideColumn>
class LocalContigMT {

public:

  LocalContigMT( const ContigSize_t contig_size):  contig_size_(contig_size), lock_strategy_(contig_size) {

    auto read_block = new (std::nothrow) NucleotideReadCount_t[(contig_size * NucleotideColumn::NUCLEOTIDE_COLUMNS)];

    data_ptr_.reset(read_block);

    if (data_ptr_.get() == nullptr) {

      ExecEnv::log().critical("Failed to allocate memory for contig. data block size: {}, nucleotides: {}",
                              contig_size, NucleotideColumn::NUCLEOTIDE_COLUMNS);

    }

    initialize(0);  // Ensure the array is properly initialized.

  }
  virtual ~LocalContigMT() = default;

  // Threadsafe read count increment.
  inline void incrementCount(const ContigOffset_t contig_offset, const Nucleotide_t nucleotide) {

    std::size_t column =  NucleotideColumn::nucleotideToColumn(nucleotide);

    if (contig_offset >= contig_size_) {

      ExecEnv::log().critical("Invalid access in incrementCount(); Contig index: {} >= Contig size: {}",
                              contig_offset, contig_size_);

    }
    if (column >= NucleotideColumn::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::NUCLEOTIDE_COLUMNS;
      ExecEnv::log().critical("incrementCount() invalid access; Nucleotide column index: {} >= Number Nucleotides: {}",
                              column, nucleotide_columns);

    }

    // pointer arithmetic, stride is NucleotideColumn::NUCLEOTIDE_COLUMNS
    NucleotideReadCount_t* access_ptr = data_ptr_.get()
                                        + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);

    lock_strategy_.incrementCount(*access_ptr, contig_offset);

  }

  inline void incrementDelete(const ContigOffset_t contig_offset) {

    incrementCount(contig_offset, NucleotideColumn::DELETE_NUCLEOTIDE);

  }

  // Reads are not locked. Do not call while updating read counts.
  inline const NucleotideReadCount_t readCount( const ContigOffset_t contig_offset,
                                                const Nucleotide_t nucleotide) const {

    std::size_t column = NucleotideColumn::nucleotideToColumn(nucleotide);

    if (contig_offset >= contig_size_) {

      ExecEnv::log().critical("readCount() invalid access; Contig index: {} >= Contig size: {}",
                              contig_offset, contig_size_);

    }
    if (column >= NucleotideColumn::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::NUCLEOTIDE_COLUMNS;
      ExecEnv::log().critical("readCount() invalid access; Nucleotide column index: {} >= Number Nucleotides: {}",
                              column, nucleotide_columns);

    }

    // pointer arithmetic, stride is NucleotideColumn::NUCLEOTIDE_COLUMNS
    const NucleotideReadCount_t  *access_ptr = data_ptr_.get()
                                               + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);
    return *access_ptr;  // read access - no mutex.

  }

  // Reads are not locked. Do not call while updating read counts.
  inline const NucleotideReadCount_t* readCount( const ContigOffset_t contig_offset) const {


    if (contig_offset >= contig_size_) {

      ExecEnv::log().critical("readCount() invalid access; Contig index: {} >= Contig size: {}",
                              contig_offset, contig_size_);

    }

    // pointer arithmetic, stride is NucleotideColumn::NUCLEOTIDE_COLUMNS
    const NucleotideReadCount_t  *access_ptr = data_ptr_.get()
                                               + (contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS);
    return access_ptr;  // read access - no mutex.

  }


  inline const ContigSize_t contigSize() const { return contig_size_; }

private:

  std::unique_ptr<NucleotideReadCount_t> data_ptr_;
  const ContigSize_t contig_size_;      // Number of of NUCLEOTIDE_COLUMNS in the contig.
  mutable std::mutex mutex_;  // Used to initialize the data.
  LockStrategy lock_strategy_; // The template multi-thread locking strategy.

  void initialize(const NucleotideReadCount_t initial_value) {

    // Initialization happens before threads are spawned but we lock the data just to make sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    for (ContigOffset_t contig_offset = 0; contig_offset < contig_size_; ++contig_offset) {
      for (ContigOffset_t column = 0; column < NucleotideColumn::NUCLEOTIDE_COLUMNS; ++column) {
        // pointer arithmetic, stride is expected_columns
        NucleotideReadCount_t *access_ptr = data_ptr_.get()
                                            + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);
        (*access_ptr) = initial_value;
      }

    }

  }

};

}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_MT_CONTIG_H
