//
// Created by kellerberrin on 19/4/20.
//

#ifndef KGL_VARIANT_FILE_IMPL_H
#define KGL_VARIANT_FILE_IMPL_H


#include <string>
#include <vector>
#include <fstream>
#include <thread>

#include "kel_exec_env.h"
#include "kel_mt_queue.h"
#include "kgl_genome_types.h"

namespace kellerberrin::genome {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plug one of the superclasses (defined in the implementation file) to read text or gzipped files.

class BaseStreamIO {

public:

  BaseStreamIO() = default;
  virtual ~BaseStreamIO() = default;

  virtual bool open(const std::string &file_name) = 0;
  virtual bool readLine(std::string &text_line) = 0;

private:

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IO object (multi-threaded) - enqueues file line reads for further processing
// The code checks the file extension, if '.vcf' then it assumes a text file, if '.gz' then it assumes a gzipped file.

class FileVCFIO  {

public:

  explicit FileVCFIO(const std::string& vcf_file_name) : vcf_file_name_(vcf_file_name),
                                                         raw_io_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  virtual ~FileVCFIO();

  void commenceIO(size_t reader_threads); // Begin reading records, spawns threads.

  [[nodiscard]] const std::string& fileName() const { return vcf_file_name_; }

  // returns an IO stream.
  [[nodiscard]] std::unique_ptr<BaseStreamIO> getSynchStream();

  // Called by each reader thread
  [[nodiscard]] std::unique_ptr<std::string> readIORecord() { return raw_io_queue_.waitAndPop(); }


private:

  std::string vcf_file_name_;
  size_t reader_threads_;
  BoundedMtQueue<std::unique_ptr<std::string>> raw_io_queue_; // The raw IO queue
  std::unique_ptr<std::thread> raw_io_thread_ptr_;

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records
  constexpr static const char* GZ_FILE_EXTENSTION_ = ".GZ"; // gzipped file assumed.

  void rawVCFIO(std::unique_ptr<BaseStreamIO>&& vcf_stream); // read/decompress from disk.

};



} // end namespace

#endif //KGL_VARIANT_FILE_IMPL_H
