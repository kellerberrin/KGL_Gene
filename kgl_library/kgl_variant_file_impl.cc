//
// Created by kellerberrin on 18/4/20.
//


#include <kel_exec_env.h>
#include "kgl_variant_file_impl.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <string>
#include <vector>
#include <memory>
#include <thread>

namespace bio = boost::iostreams;
// Implementation file classes need to be defined within namespaces.
namespace kellerberrin::genome {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Text IO.
class TextStreamIO : public BaseStreamIO {

public:

  TextStreamIO() = default;

  TextStreamIO(const TextStreamIO &) = delete;

  ~TextStreamIO() override = default;

  bool open(const std::string &file_name) override;

  bool readLine(std::string &text_line) override { return not std::getline(file_, text_line).eof(); }

private:

  std::ifstream file_;

};


bool TextStreamIO::open(const std::string &file_name) {

  try {

    // Open input file.

    file_.open(file_name);
    if (not file_.good()) {

      ExecEnv::log().error("TextStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }
  }
  catch (std::exception const &e) {

    ExecEnv::log().error("TextStreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gzip IO

class GZStreamIO : public BaseStreamIO {

public:

  GZStreamIO() = default;

  GZStreamIO(const GZStreamIO &) = delete;

  ~GZStreamIO() override = default;

  bool open(const std::string &file_name) override;

  bool readLine(std::string &text_line) override { return not std::getline(gz_file_, text_line).eof(); }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream gz_file_;

};

bool GZStreamIO::open(const std::string &file_name) {

  try {

    // Open input file.

    file_.open(file_name, std::ios_base::in | std::ios_base::binary);

    if (not file_.good()) {

      ExecEnv::log().error("GZStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }

    gz_file_.push(bio::gzip_decompressor());
    gz_file_.push(file_);

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("GZStreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

FileVCFIO::~FileVCFIO() {

  if (raw_io_thread_ptr_) raw_io_thread_ptr_->join();

}

std::unique_ptr<BaseStreamIO> FileVCFIO::getSynchStream() {

  std::string file_ext = Utility::toupper(Utility::fileExtension(fileName()));

  if (file_ext == GZ_FILE_EXTENSTION_) {

    return std::make_unique<GZStreamIO>();

  } else {   // Assume a text file.


    return std::make_unique<TextStreamIO>();

  }

}

void FileVCFIO::commenceIO(size_t reader_threads) {

  reader_threads_ = reader_threads; // the number of consumers; used for shutdown.
  raw_io_thread_ptr_ = std::make_unique<std::thread>(&FileVCFIO::rawVCFIO, this, getSynchStream());

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pass in the text or gzip stream as an rvalue to a pointer.
void FileVCFIO::rawVCFIO(std::unique_ptr<BaseStreamIO> &&vcf_stream) {

  try {

    if (not vcf_stream->open(fileName())) {

      ExecEnv::log().critical("FileVCFIO; I/O error; could not open VCF file: {}", fileName());

    }

    while (true) {

      std::unique_ptr<std::string> line_record_ptr = std::make_unique<std::string>();

      if (not vcf_stream->readLine(*line_record_ptr)) {

        // Enqueue the null eof indicator for each consumer thread.
        for (size_t i = 0; i < reader_threads_; ++i) {

          raw_io_queue_.push(std::unique_ptr<std::string>(nullptr));

        }

        // Terminate the read line loop.
        break;

      }

      // Check we have a valid line record.
      if (line_record_ptr->empty()) {

        ExecEnv::log().error("FileVCFIO; Empty VCF record returned");
        continue;

      }

      raw_io_queue_.push(std::move(line_record_ptr));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("FileVCFIO; VCF file unexpected I/O exception: {}", e.what());

  }

}



} // namespace