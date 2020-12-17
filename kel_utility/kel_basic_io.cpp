//
// Created by kellerberrin on 11/5/20.
//

#include <kel_exec_env.h>
#include "kel_basic_io.h"
#include "kel_utility.h"
#include "kel_bzip.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <string>
#include <vector>
#include <memory>

namespace bio = boost::iostreams;
// Implementation file classes need to be defined within namespaces.
namespace kellerberrin {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Text IO.
class TextStreamIO : public BaseStreamIO {

public:

  TextStreamIO() = default;

  TextStreamIO(const TextStreamIO &) = delete;

  ~TextStreamIO() override = default;

  bool open(const std::string &file_name) override;

  IOLineRecord readLine() override {

    std::unique_ptr<std::string> line_text_ptr(std::make_unique<std::string>());
    if (not std::getline(file_, *line_text_ptr).eof()) {

      ++record_counter_;
      return IOLineRecord(std::pair<size_t, std::unique_ptr<std::string>>(record_counter_, std::move(line_text_ptr)));

    } else {

      return std::nullopt;

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

    std::unique_ptr<std::string> line_text_ptr(std::make_unique<std::string>());
    if (not std::getline(gz_file_, *line_text_ptr).eof()) {

      ++record_counter_;
      return IOLineRecord(std::pair<size_t, std::unique_ptr<std::string>>(record_counter_, std::move(line_text_ptr)));

    } else {

      return std::nullopt;

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


std::optional<std::unique_ptr<BaseStreamIO>> BaseStreamIO::getReaderStream(const std::string &file_name) {

  std::string file_ext = Utility::toupper(Utility::fileExtension(file_name));

  // Open a compressed gzipped file based on the file extension.
  if (file_ext == GZ_FILE_EXTENSTION_)  {


    std::unique_ptr<BaseStreamIO> gz_stream(std::make_unique<GZStreamIO>());
    if (not gz_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return gz_stream;

    }

  } else if (file_ext == BGZ_FILE_EXTENSTION_) {

    std::unique_ptr<BaseStreamIO> gz_stream(std::make_unique<GZBlockDecompression>());
    if (not gz_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return gz_stream;

    }

  } else {   // Assume a text file.

    std::unique_ptr<BaseStreamIO> text_stream(std::make_unique<TextStreamIO>());
    if (not text_stream->open(file_name)) {

      return std::nullopt;

    } else {

      return text_stream;

    }

  }

}



} // namespace