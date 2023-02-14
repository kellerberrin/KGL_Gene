//
// Created by kellerberrin on 14/02/23.
//

#ifndef KEL_FILE_IO_H
#define KEL_FILE_IO_H

#include "kel_basic_io.h"

#include <string>
#include <memory>
#include <fstream>


namespace kellerberrin {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plain text IO.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TextStreamIO : public BaseStreamIO {

public:

  TextStreamIO() = default;
  TextStreamIO(const TextStreamIO &) = delete;
  ~TextStreamIO() override { close(); };

  [[nodiscard]] bool open(const std::string &file_name) override;
  [[nodiscard]] IOLineRecord readLine() override;
  void close() override { record_counter_ = 0; file_.close(); }

  [[nodiscard]] static std::optional<std::unique_ptr<BaseStreamIO>> getStreamIO( const std::string& file_name);

private:

  std::ifstream file_;
  size_t record_counter_{0};

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of gzip '.gz' decompression uses Pimpl idiom.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GZStreamIOImpl; // Currently implemented using boost::streams.
class GZStreamIO : public BaseStreamIO {

public:

  GZStreamIO();
  GZStreamIO(const GZStreamIO &) = delete;
  ~GZStreamIO() override;

  // Pimpl redirected functions.
  [[nodiscard]] bool open(const std::string &file_name) override;
  [[nodiscard]] IOLineRecord readLine() override;
  void close() override;

  [[nodiscard]] static std::optional<std::unique_ptr<BaseStreamIO>> getStreamIO( const std::string& file_name);

private:

  std::unique_ptr<GZStreamIOImpl> pimpl_streamio_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of Burrows-Wheeler '.bz2' decompression uses Pimpl idiom.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BZ2StreamIOImpl; // Currently implemented using boost::streams.
class BZ2StreamIO : public BaseStreamIO {

public:

  BZ2StreamIO();
  BZ2StreamIO(const BZ2StreamIO &) = delete;
  ~BZ2StreamIO() override;

  [[nodiscard]] bool open(const std::string &file_name) override;
  [[nodiscard]] IOLineRecord readLine() override;
  void close() override;

  [[nodiscard]] static std::optional<std::unique_ptr<BaseStreamIO>> getStreamIO( const std::string& file_name);

private:

  std::unique_ptr<BZ2StreamIOImpl> pimpl_streamio_;

};



} // namespace


#endif //KEL_FILE_IO_H
