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

  NumpyContigMT(  NucleotideReadCount_t *data_ptr,
                  const ContigSize_t contig_size,
                  const std::size_t nucleotides) :  data_ptr_(data_ptr),
                                                    contig_size_(contig_size),
                                                    lock_strategy_(contig_size) {

    if (nucleotides != NucleotideColumn::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = NucleotideColumn::NUCLEOTIDE_COLUMNS;
      ExecEnv::log().critical("Numpy array expected to have: {} NUCLEOTIDE_COLUMNS (columns), columns specified: {}",
                              nucleotide_columns, nucleotides);

    }

    initialize(0);  // Ensure the array is properly initialized.

  }
  virtual ~NumpyContigMT() = default;

  // Threadsafe read count increment.
  inline void incrementCount(const ContigOffset_t contig_offset, std::size_t column) {

    if (contig_offset >= contig_size_) {

      ExecEnv::log().critical("Invalid access in incrementCount(); Contig index: {} >= Contig size: {}",
                              contig_offset, contig_size_);

    }
    if (column >= ExtendDNA5::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = ExtendDNA5::NUCLEOTIDE_COLUMNS;
      ExecEnv::log().critical("incrementCount() invalis access; Nucleotide column index: {} >= Number Nucleotides: {}",
                              column, nucleotide_columns);

    }

    // pointer arithmetic, stride is NucleotideColumn::NUCLEOTIDE_COLUMNS
    NucleotideReadCount_t* access_ptr = data_ptr_ + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);

    lock_strategy_.incrementCount(*access_ptr, contig_offset);

  }

  // Reads are not locked. Do not call while updating read counts.
  inline const NucleotideReadCount_t readCount( const ContigOffset_t contig_offset,
                                                std::size_t column) const {

    if (contig_offset >= contig_size_) {

      ExecEnv::log().critical("Invalid access in incrementCount(); Contig index: {} >= Contig size: {}",
                              contig_offset, contig_size_);

    }
    if (column >= ExtendDNA5::NUCLEOTIDE_COLUMNS) {

      ContigOffset_t nucleotide_columns = ExtendDNA5::NUCLEOTIDE_COLUMNS;
      ExecEnv::log().critical("readCount() invalid access; Nucleotide column index: {} >= Number Nucleotides: {}",
                              column, nucleotide_columns);

    }

    // pointer arithmetic, stride is NucleotideColumn::NUCLEOTIDE_COLUMNS
    const NucleotideReadCount_t  *access_ptr = data_ptr_
                                               + ((contig_offset * NucleotideColumn::NUCLEOTIDE_COLUMNS) + column);
    return *access_ptr;  // read access - no mutex.

  }

  inline const ContigSize_t contigSize() const { return contig_size_; }

private:

  // Yes folks, that's a raw pointer (watch it snort and shake).
  NucleotideReadCount_t *data_ptr_; // Snort, shake, snort, shake ...
  const ContigSize_t contig_size_;      // Number of of NUCLEOTIDE_COLUMNS in the contig.
  mutable std::mutex mutex_;  // Used to initialize the data.
  LockStrategy lock_strategy_; // The template multi-thread locking strategy.

  void initialize(const NucleotideReadCount_t initial_value) {

    // Initialization happens before threads are spawned but we lock the data just to make sure.
    std::lock_guard<std::mutex> lock(mutex_) ;

    for (ContigOffset_t contig_offset = 0; contig_offset < contig_size_; ++contig_offset) {
      for (ContigOffset_t column = 0; column < ExtendDNA5::NUCLEOTIDE_COLUMNS; ++column) {
        // pointer arithmetic, stride is expected_columns
        NucleotideReadCount_t *access_ptr = data_ptr_
                                            + ((contig_offset * ExtendDNA5::NUCLEOTIDE_COLUMNS) + column);
        (*access_ptr) = initial_value;
      }

    }

  }

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_MT_NUMPY_H
