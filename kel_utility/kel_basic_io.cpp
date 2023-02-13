//
// Created by kellerberrin on 11/5/20.
//

#include <kel_exec_env.h>
#include "kel_basic_io.h"
#include "kel_utility.h"
#include "kel_bzip_workflow.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <string>
#include <vector>
#include <memory>
#include <fstream>


namespace bio = boost::iostreams;
// Implementation file classes need to be defined within namespaces.
namespace kellerberrin {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// After moving data out of the object, it is empty() but not EOF().
std::pair<size_t, std::string> IOLineRecord::getLineData() {

  if (empty() or EOFRecord()) {

    // Complain and return the empty string.
    ExecEnv::log().warn("IOLineRecord::getLineData; attempt to move data from an empty()/EOFRecord() object");
    return {0, {}};

  }

  size_t out_line = line_count_;
  line_count_ = 0;
  empty_ = true;
  return {out_line, std::move(line_data_)};

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Text IO.
class TextStreamIO : public BaseStreamIO {

public:

  TextStreamIO() = default;

  TextStreamIO(const TextStreamIO &) = delete;

  ~TextStreamIO() override = default;

  bool open(const std::string &file_name) override;

  IOLineRecord readLine() override {

    std::string line_text;
    if (not std::getline(file_, line_text).eof()) {

      ++record_counter_;
      return IOLineRecord(record_counter_, std::move(line_text));

    } else {

      return IOLineRecord::createEOFMarker();

    }

  }

private:

  std::ifstream file_;

};


bool TextStreamIO::open(const std::string &file_name) {

  try {

    record_counter_ = 0;
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

  inline IOLineRecord readLine() override {

    std::string line_text;
    if (not std::getline(gz_file_, line_text).eof()) {

      ++record_counter_;
      return IOLineRecord(record_counter_, std::move(line_text));

    } else {

      return IOLineRecord::createEOFMarker();

    }

  }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream gz_file_;

};

bool GZStreamIO::open(const std::string &file_name) {

  try {

    record_counter_ = 0;
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gzip IO

class BZ2StreamIO : public BaseStreamIO {

public:

  BZ2StreamIO() = default;

  BZ2StreamIO(const BZ2StreamIO &) = delete;

  ~BZ2StreamIO() override = default;

  bool open(const std::string &file_name) override;

  inline IOLineRecord readLine() override {

    std::string line_text;
    if (not std::getline(bz2_file_, line_text).eof()) {

      ++record_counter_;
      return IOLineRecord(record_counter_, std::move(line_text));

    } else {

      return IOLineRecord::createEOFMarker();

    }

  }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream bz2_file_;

};

bool BZ2StreamIO::open(const std::string &file_name) {

  try {

    record_counter_ = 0;
    // Open input file.

    file_.open(file_name, std::ios_base::in | std::ios_base::binary);

    if (not file_.good()) {

      ExecEnv::log().error("GZStreamIO; I/O error; could not open file: {}", file_name);
      return false;

    }

    bz2_file_.push(bio::bzip2_decompressor());
    bz2_file_.push(file_);

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("BZ2StreamIO; Opening file: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  return true;

}


std::optional<std::unique_ptr<BaseStreamIO>> BaseStreamIO::getStreamIO(const std::string &file_name, size_t decompression_threads) {

  std::string file_ext = Utility::toupper(Utility::fileExtension(file_name));

  // Open a compressed gzipped file based on the file extension.
  if (file_ext == GZ_FILE_EXTENSTION_) {

    std::unique_ptr<BaseStreamIO> gz_stream;
    if (BGZStream::verify(file_name)) { // Check if block gzipped.

      ExecEnv::log().info("File structure verified as bgz format, parser uses bgz reader.");
      gz_stream = std::make_unique<BGZStream>(decompression_threads);

    } else {

      ExecEnv::log().info("File structure is not in bgz format, parser uses a general purpose gzip reader.");
      gz_stream = std::make_unique<GZStreamIO>();

    }
    if (not gz_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return gz_stream;

    }

  } else if (file_ext == BGZ_FILE_EXTENSTION_) { // .bgz extension is just assumed to be block gzipped.

    std::unique_ptr<BaseStreamIO> gz_stream(std::make_unique<BGZStream>(decompression_threads));
    if (not gz_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return gz_stream;

    }

  } else if (file_ext == BZ2_FILE_EXTENSTION_) { // If .bz2 then use Burrows Wheeler decompression.

    std::unique_ptr<BaseStreamIO> bz2_stream(std::make_unique<BZ2StreamIO>());

    if (not bz2_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return bz2_stream;

    }

  } else {   // Otherwise assume a text file.

    std::unique_ptr<BaseStreamIO> text_stream(std::make_unique<TextStreamIO>());
    if (not text_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return text_stream;

    }

  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

// Uses the filename extension heuristic documented above to open the underlying StreamIO
// The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
bool MTStreamIO::open(const std::string &file_name, size_t decompression_threads) {

  auto stream_opt = BaseStreamIO::getStreamIO(file_name, decompression_threads);
  if (not stream_opt) {

    ExecEnv::log().error("MTStreamIO::open; could not open stream for file: {}", file_name);
    return false;

  }

  EOF_received_ = false;
  stream_ptr_ = std::move(stream_opt.value());
  line_io_thread_.queueThreads(WORKER_THREAD_COUNT);
  line_io_thread_.enqueueVoid(&MTStreamIO::enqueueIOLineRecord, this);

  return true;

}

void MTStreamIO::close() {

  EOF_received_ = true;
  line_io_thread_.joinThreads();
  stream_ptr_ = nullptr;
  line_io_queue_.clear();

}

void MTStreamIO::enqueueIOLineRecord() {

  while(true) {

    auto line_record = stream_ptr_->readLine();

    if (line_record.EOFRecord()) {

      EOF_received_ = true;
      line_io_queue_.push(std::move(line_record));
      break;

    }

    line_io_queue_.push(std::move(line_record));

  }

}

IOLineRecord MTStreamIO::readLine() {

  // Dont block on EOF
  std::lock_guard<std::mutex> lock(mutex_);
  if (EOF_received_) return IOLineRecord::createEOFMarker();
  return line_io_queue_.waitAndPop();

}



} // namespace
