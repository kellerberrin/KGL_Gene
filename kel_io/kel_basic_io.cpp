//
// Created by kellerberrin on 11/5/20.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kel_basic_io.h"
#include "kel_file_io.h"
#include "kel_bzip_workflow.h"


// Implementation file classes need to be defined within namespaces.
namespace kellerberrin {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// After moving data out of the object, it is empty() but not EOF().
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns an IO stream that is either a normal stream or a compressed stream based on the file name extension.
// If '.bgz' then it uses a multi-threaded and memory efficient algorithm to decompress as it reads.
// If '.gz' then the file is analyzed to see if it is in '.bgz' format.
// If not in '.bgz' format then a standard (boost) single thread decompression algorithm is used.
// If '.bz2' then Burrows-Wheeler compression is assumed.
// If the file is not one of the above types, it is assumed to be an uncompressed record based text file.
// Note the stream is returned open and ready for processing, std::nullopt is returned if there is a problem.
// The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::optional<std::unique_ptr<BaseStreamIO>> BaseStreamIO::getStreamIO(const std::string &file_name) {

  std::string file_ext = Utility::toupper(Utility::fileExtension(file_name));

  // Open a compressed gzipped file based on the file extension.
  if (file_ext == GZ_FILE_EXTENSTION_) {

    std::optional<std::unique_ptr<BaseStreamIO>> gz_stream_opt;
    if (BGZStreamIO::verify(file_name)) { // Check if block gzipped.

      ExecEnv::log().info("File structure verified as bgz format, parser uses bgz reader.");
      return BGZStreamIO::getStreamIO(file_name);

    } else {

      ExecEnv::log().info("File structure is not in bgz format, parser uses a general purpose gzip reader.");
      return GZStreamIO::getStreamIO(file_name);

    }

  } else if (file_ext == BGZ_FILE_EXTENSTION_) { // .bgz extension is just assumed to be block gzipped.

    return BGZStreamIO::getStreamIO(file_name);

  } else if (file_ext == BZ2_FILE_EXTENSTION_) { // If .bz2 then use Burrows Wheeler decompression.

    return BZ2StreamIO::getStreamIO(file_name);

  }

  // Otherwise assume a text file.
  return TextStreamIO::getStreamIO(file_name);

}


} // namespace
