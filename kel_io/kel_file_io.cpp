//
// Created by kellerberrin on 14/02/23.
//

#include "kel_file_io.h"
#include "kel_exec_env.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>


namespace bio = boost::iostreams;

namespace kellerberrin {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


IOLineRecord TextStreamIO::readLine() {

  std::string line_text;
  try {

    if (std::getline(file_, line_text)) {

      ++record_counter_;
      return IOLineRecord(record_counter_, std::move(line_text));

    } else {

      return IOLineRecord::createEOFMarker();

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("TextStreamIO::readLine; I/O exception reading file: {}", e.what());
    return IOLineRecord::createEOFMarker();

  }

}


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



std::optional<std::unique_ptr<BaseStreamIO>> TextStreamIO::getStreamIO( const std::string& file_name) {

  auto stream_ptr = std::make_unique<TextStreamIO>();
  if (stream_ptr->open(file_name)) {

    return stream_ptr;

  }

  ExecEnv::log().error("TextStreamIO::getStreamIO; error opening file: {}", file_name);
  return std::nullopt;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of gzip '.gz' and bzip2 '.bz2' decompression uses boost::iostreams::filtering_istream.
// The two stream types differ only in the boost filter, so a single template implementation is shared.
// GZStreamIOImpl/BZ2StreamIOImpl remain opaque forward-declared classes in the header (no boost exposure).
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename Filter>
class FilterStreamIOImpl {

public:

  FilterStreamIOImpl() = default;
  FilterStreamIOImpl(const FilterStreamIOImpl &) = delete;
  ~FilterStreamIOImpl() = default;

  bool open(const std::string &file_name, const char* stream_name) {

    try {

      record_counter_ = 0;
      // Open input file.

      file_.open(file_name, std::ios_base::in | std::ios_base::binary);

      if (not file_.good()) {

        ExecEnv::log().error("{}; I/O error; could not open file: {}", stream_name, file_name);
        return false;

      }

      filter_file_.push(Filter());
      filter_file_.push(file_);

    }
    catch (std::exception const &e) {

      ExecEnv::log().error("{}; Opening file: {} unexpected I/O exception: {}", stream_name, file_name, e.what());
      return false;

    }

    return true;

  }

  void close() { record_counter_ = 0; filter_file_.reset(); }

  IOLineRecord readLine() {

    std::string line_text;
    try {

      if (std::getline(filter_file_, line_text)) {

        ++record_counter_;
        return IOLineRecord(record_counter_, std::move(line_text));

      } else {

        return IOLineRecord::createEOFMarker();

      }

    }
    catch (std::exception const &e) {

      ExecEnv::log().error("FilterStreamIOImpl::readLine; I/O exception reading file: {}", e.what());
      return IOLineRecord::createEOFMarker();

    }

  }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream filter_file_;
  size_t record_counter_{0};

};


// Concrete Pimpl types, matching the opaque forward declarations in the header.
class GZStreamIOImpl final : public FilterStreamIOImpl<bio::gzip_decompressor> {};
class BZ2StreamIOImpl final : public FilterStreamIOImpl<bio::bzip2_decompressor> {};


// Pimpl redirection.
GZStreamIO::GZStreamIO() { pimpl_streamio_ = std::make_unique<GZStreamIOImpl>(); }
GZStreamIO::~GZStreamIO() { close(); }
bool GZStreamIO::open(const std::string &file_name) { return pimpl_streamio_->open(file_name, "GZStreamIO"); }
IOLineRecord GZStreamIO::readLine() { return pimpl_streamio_->readLine(); }
void GZStreamIO::close() { pimpl_streamio_->close(); }


std::optional<std::unique_ptr<BaseStreamIO>> GZStreamIO::getStreamIO( const std::string& file_name) {

  auto stream_ptr = std::make_unique<GZStreamIO>();
  if (stream_ptr->open(file_name)) {

    return stream_ptr;

  }

  ExecEnv::log().error("GZStreamIO::getStreamIO; error opening file: {}", file_name);
  return std::nullopt;

}


// Pimpl redirection.
BZ2StreamIO::BZ2StreamIO() { pimpl_streamio_ = std::make_unique<BZ2StreamIOImpl>(); }
BZ2StreamIO::~BZ2StreamIO() { close(); }
bool BZ2StreamIO::open(const std::string &file_name) { return pimpl_streamio_->open(file_name, "BZ2StreamIO"); }
IOLineRecord BZ2StreamIO::readLine() { return pimpl_streamio_->readLine(); }
void BZ2StreamIO::close() { pimpl_streamio_->close(); }


std::optional<std::unique_ptr<BaseStreamIO>> BZ2StreamIO::getStreamIO( const std::string& file_name) {

  auto stream_ptr = std::make_unique<BZ2StreamIO>();
  if (stream_ptr->open(file_name)) {

    return stream_ptr;

  }

  ExecEnv::log().error("BZ2StreamIO::getStreamIO; error opening file: {}", file_name);
  return std::nullopt;

}

} // namespace
