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
  if (not std::getline(file_, line_text).eof()) {

  ++record_counter_;
  return IOLineRecord(record_counter_, std::move(line_text));

  } else {

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
// Implementation of gzip '.gz' decompression uses uses boost::iostreams::filtering_istream.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GZStreamIOImpl {

public:

  GZStreamIOImpl() = default;

  GZStreamIOImpl(const GZStreamIOImpl &) = delete;

  ~GZStreamIOImpl() = default;


  bool open(const std::string &file_name);

  void close() { record_counter_ = 0; gz_file_.reset(); }

  IOLineRecord readLine();

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream gz_file_;
  size_t record_counter_{0};

};


IOLineRecord GZStreamIOImpl::readLine() {

  std::string line_text;
  if (not std::getline(gz_file_, line_text).eof()) {

    ++record_counter_;
    return IOLineRecord(record_counter_, std::move(line_text));

  } else {

    return IOLineRecord::createEOFMarker();

  }

}


bool GZStreamIOImpl::open(const std::string &file_name) {

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

// Pimpl redirection.
GZStreamIO::GZStreamIO() { pimpl_streamio_ = std::make_unique<GZStreamIOImpl>(); }
GZStreamIO::~GZStreamIO() { close(); }
bool GZStreamIO::open(const std::string &file_name) { return pimpl_streamio_->open(file_name); }
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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of Burrows-Wheeler '.bz2' decompression uses boost::iostreams::filtering_istream.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BZ2StreamIOImpl  {

public:

  BZ2StreamIOImpl() = default;
  BZ2StreamIOImpl(const BZ2StreamIOImpl &) = delete;
  ~BZ2StreamIOImpl() = default;

  bool open(const std::string &file_name);

  IOLineRecord readLine();

  void close() { record_counter_ = 0; bz2_file_.reset(); }

private:

  std::ifstream file_;
  boost::iostreams::filtering_istream bz2_file_;
  size_t record_counter_{0};

};


IOLineRecord BZ2StreamIOImpl::readLine() {

  std::string line_text;
  if (not std::getline(bz2_file_, line_text).eof()) {

    ++record_counter_;
    return IOLineRecord(record_counter_, std::move(line_text));

  } else {

    return IOLineRecord::createEOFMarker();

  }

}


bool BZ2StreamIOImpl::open(const std::string &file_name) {

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

// Pimpl redirection.
BZ2StreamIO::BZ2StreamIO() { pimpl_streamio_ = std::make_unique<BZ2StreamIOImpl>(); }
BZ2StreamIO::~BZ2StreamIO() { close(); }
bool BZ2StreamIO::open(const std::string &file_name) { return pimpl_streamio_->open(file_name); }
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
