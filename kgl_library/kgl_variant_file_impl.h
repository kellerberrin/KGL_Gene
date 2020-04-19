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
#include "kgl_variant_file.h"

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual IO object (multi-threaded)


class BaseVCFIO {

public:

  BaseVCFIO(const BaseVCFIO&) = delete;
  explicit BaseVCFIO(std::string vcf_file_name) : vcf_file_name_(vcf_file_name) {}
  virtual ~BaseVCFIO() = default;


  virtual void commenceIO() = 0;

  [[nodiscard]] virtual std::unique_ptr<VcfRecord> readVCFRecord() = 0;

  [[nodiscard]] virtual VcfHeaderInfo VCFReadHeader() = 0;

  [[nodiscard]] virtual const std::vector<std::string>& getGenomeNames() const = 0;

  [[nodiscard]] const std::string& fileName() const { return vcf_file_name_; }

private:

  std::string vcf_file_name_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Seqan PIMPL facade IO object (multi-threaded)

class SeqanVCFIO : public BaseVCFIO {

public:

  explicit SeqanVCFIO(const std::string& vcf_file_name);
  ~SeqanVCFIO() override;


  void commenceIO() override; // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord() override;

  [[nodiscard]] VcfHeaderInfo VCFReadHeader() override;

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const override { return parseheader_.getGenomes(); }

private:

  class SeqanVCFIOImpl;       // Forward declaration of seqan vcf read implementation class
  std::unique_ptr<SeqanVCFIOImpl> seqan_vcf_impl_ptr_;    // Seqan VCF parser PIMPL
  VCFParseHeader parseheader_;    // Get genome and contig information.

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IO object (multi-threaded) - uses standard File IO functions to parse VCF file.
// The code checks the file extension, if '.vcf' then it assumes a text file, if '.gz' then it assumes a gzipped file.

class FileVCFIO : public BaseVCFIO {

public:

  explicit FileVCFIO(const std::string& vcf_file_name) : BaseVCFIO(vcf_file_name),
                                                         raw_io_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                         vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  ~FileVCFIO() override;


  void commenceIO() override; // Begin reading records, spawns threads.

  [[nodiscard]] std::unique_ptr<VcfRecord> readVCFRecord() override;

  [[nodiscard]] VcfHeaderInfo VCFReadHeader() override;

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const override { return parseheader_.getGenomes(); }

private:

  BoundedMtQueue<std::unique_ptr<std::string>> raw_io_queue_; // The raw IO queue
  BoundedMtQueue<std::unique_ptr<VcfRecord>> vcf_record_queue_; // Parsed VCF record queue

  std::unique_ptr<std::thread> raw_io_thread_ptr_;
  std::vector<std::thread> vcf_record_thread_vec_;

  VCFParseHeader parseheader_;    // Get genome and contig information.

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records
  static constexpr const long PARSER_THREADS_{4};         // Threads parsing into vcf_records.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '#' then a header record.
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are format and genotype fields (header specified).
  static constexpr const char* VCF_FIELD_DELIMITER_{"\t"};   // VCF Field separator.
  constexpr static const char* VCF_FILE_EXTENSTION_ = ".VCF"; // Valid file extensions (case insensitive)
  constexpr static const char* GZ_FILE_EXTENSTION_ = ".GZ"; // gzipped VCF file assumed.

  template<class IOStream> void rawVCFIO(); // read/decompress from disk.
  void enqueueVCFRecord(); // enqueue vcf_records.
  bool parseVCFRecord(const std::unique_ptr<std::string>& line_record_ptr, const std::unique_ptr<VcfRecord>& vcf_record_ptr);
  bool moveToVcfRecord(std::vector<std::string>& fields, VcfRecord& vcf_record);

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as a template argument.
template<class IOStream> void FileVCFIO::rawVCFIO() {

  try {

    IOStream vcf_stream;

    if (not vcf_stream.open(fileName())) {

      ExecEnv::log().critical("I/O error; could not open VCF file: {}", fileName());

    }

    while (true) {

      std::unique_ptr<std::string> line_record_ptr = std::make_unique<std::string>();

      if (not vcf_stream.readLine(*line_record_ptr)) {

        // Enqueue the null eof indicator for each consumer thread.
        for (size_t i = 0; i < vcf_record_thread_vec_.size(); ++i) {

          raw_io_queue_.push(std::unique_ptr<std::string>(nullptr));

        }
        break;

      }

      // Check we have read a non-empty string.
      if (line_record_ptr->empty()) {

        ExecEnv::log().error("Empty VCF record string returned");
        continue;

      }

      // Skip header records.
      if ((*line_record_ptr)[0] == HEADER_CHAR_) {

        continue;

      }

      raw_io_queue_.push(std::move(line_record_ptr));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected Seqan I/O exception: {}", e.what());

  }

}


} // end namespace

#endif //KGL_VARIANT_FILE_IMPL_H
