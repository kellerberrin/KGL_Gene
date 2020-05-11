//
// Created by kellerberrin on 19/4/20.
//

#ifndef KGL_VARIANT_FILE_IMPL_H
#define KGL_VARIANT_FILE_IMPL_H



#include "kel_exec_env.h"
#include "kel_mt_queue.h"
#include "kel_basic_io.h"
#include "kgl_genome_types.h"

#include <string>
#include <vector>
#include <fstream>
#include <thread>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IO object (multi-threaded) - enqueues file line reads for further processing

class FileVCFIO  {

public:

  explicit FileVCFIO(const std::string& vcf_file_name) : vcf_file_name_(vcf_file_name),
                                                         raw_io_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  virtual ~FileVCFIO();

  // Begin reading records, spawns thread(s). Specifies how many reader threads.
  void commenceIO(size_t reader_threads);

  [[nodiscard]] const std::string& fileName() const { return vcf_file_name_; }

  // Called by each reader thread
  [[nodiscard]] IOLineRecord readIORecord() { return raw_io_queue_.waitAndPop(); }


private:

  std::string vcf_file_name_;
  size_t reader_threads_;
  BoundedMtQueue<IOLineRecord> raw_io_queue_; // The raw tidal IO queue
  std::unique_ptr<std::thread> raw_io_thread_ptr_;

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records

  void rawVCFIO(); // The read/decompress from disk and place on the tidal queue.

};



} // end namespace

#endif //KGL_VARIANT_FILE_IMPL_H
