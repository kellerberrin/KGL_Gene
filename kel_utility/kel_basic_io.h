//
// Created by kellerberrin on 11/5/20.
//

#ifndef KEL_BASIC_IO_H
#define KEL_BASIC_IO_H

#include <string>
#include <vector>
#include <optional>
#include <fstream>

#include "kel_exec_env.h"

namespace kellerberrin {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plug one of the superclasses (defined in the implementation file) to read text or gzipped files.

// Returns an optional pair of line count and a pointer to the line data.
// The data is absent on EOF or read error.
using IOLineRecord = std::optional<std::pair<size_t, std::unique_ptr<std::string>>>;
// The eof marker.
constexpr static const std::nullopt_t QUEUED_EOF_MARKER = std::nullopt;


class BaseStreamIO {

public:

  BaseStreamIO() = default;
  virtual ~BaseStreamIO() = default;

  virtual bool open(const std::string &file_name) = 0;
  virtual IOLineRecord readLine() = 0;

  // Returns an IO stream that is either a normal stream or a compressed stream based on the file name extension.
  // If '.gz' or '.bgz' then it assumes a gzipped file and decompresses as it reads.
  // Note the stream is returned open and ready for processing, std::nullopt is returned if there is a problem.
  [[nodiscard]] static std::optional<std::unique_ptr<BaseStreamIO>> getReaderStream(const std::string& file_name);

protected:

  size_t record_counter_{0};

  constexpr static const char* GZ_FILE_EXTENSTION_ = ".GZ"; // gzipped file assumed.
  constexpr static const char* BGZ_FILE_EXTENSTION_ = ".BGZ"; // gzipped file assumed.
  constexpr static const char* BZ2_FILE_EXTENSTION_ = ".BZ2"; // Burrows-Wheeler compression assumed.

};


} // namespace


#endif //KEL_BASIC_IO_H
