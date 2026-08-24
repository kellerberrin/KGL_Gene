//
// Created by kellerberrin on 11/5/20.
//

#ifndef KEL_BASIC_IO_H
#define KEL_BASIC_IO_H

#include <string>
#include <memory>
#include <string_view>
#include <optional>
#include <utility>


namespace kellerberrin {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A lightweight object that contains the readLine() output of StreamIO with file line number and associated line text.
// Note that the line text has been stripped of the terminating '\n'.
// Once the line data has been moved from the object using getLineData() it is tagged as empty().
// A subsequent attempt to access the moved data using getLineData() will generate a warning.
// A static function createEOFMarker() creates an object to serve as EOF, this can also be pushed onto a queue.
// Since the stored std::string data can be large, for optimal performance the object cannot be copied.
// The object is designed to be pushed and popped from thread safe queues without incurring significant processing overhead.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// A lightweight record containing a line number and associated line text from readLine().
/// The line text is stripped of the terminating '\n'. After getLineData() moves the data out, the object is empty().
/// Use createEOFMarker() to create an EOF sentinel. The object cannot be copied but can be moved.
class IOLineRecord {

public:

  IOLineRecord(size_t line_count, std::string&& line_data) noexcept : line_count_(line_count), line_data_(std::move(line_data)) {}

  IOLineRecord(const IOLineRecord&) = delete;
  IOLineRecord(IOLineRecord&&) noexcept = default;
  IOLineRecord& operator=(const IOLineRecord&) = delete;
  IOLineRecord& operator=(IOLineRecord&&) noexcept = default;
  ~IOLineRecord() = default;

  /// Moves the line data out of the object. After this call the object is empty() but not EOF().
  /// A subsequent call on an empty() or EOFRecord() object logs a warning and returns {0, ""}.
  [[nodiscard]] std::pair<size_t, std::string> getLineData();
  [[nodiscard]] size_t lineCount() const noexcept { return line_count_; }
  [[nodiscard]] std::string_view getView() const noexcept { return {line_data_}; }
  [[nodiscard]] bool EOFRecord() const noexcept { return EOF_; }
  [[nodiscard]] bool empty() const noexcept { return empty_; }

  /// Returns an EOF marker suitable for pushing onto a queue.
  [[nodiscard]] static IOLineRecord createEOFMarker() { return {}; }

private:

  IOLineRecord() : EOF_{true} {} // Only create as an EOF marker.

  size_t line_count_{0}; // Actual line counts begin at 1.
  std::string line_data_; // This string is NOT '\n' terminated.
  bool EOF_{false};
  bool empty_{false};

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plug one of the superclasses (defined in the implementation file) to read text or gzipped files.
// All decompression is done transparently behind this object.
// In particular large '.bgz' are decompressed using multiple threads, this is useful for large VCF files.
// Note that the readline() virtual function is not, in general, multithreaded.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Abstract base class for stream-based line readers.
/// Use getStreamIO() to obtain a concrete implementation based on the file name extension.
/// Decompression (bgz, gz, bz2) is handled transparently by the returned implementation.
class BaseStreamIO {

public:

  BaseStreamIO() = default;
  virtual ~BaseStreamIO() = default;

  virtual bool open(const std::string &file_name) = 0;
  virtual void close() = 0;
  virtual IOLineRecord readLine() = 0;

  /// Returns an IO stream appropriate for the file name extension.
  /// '.bgz' uses multi-threaded decompression; '.gz' is checked for bgz format then falls back to standard gzip;
  /// '.bz2' uses Burrows-Wheeler decompression; otherwise the file is treated as uncompressed text.
  /// The stream is returned open and ready for processing. Returns std::nullopt on failure.
  /// The decompression_threads argument is only used by the bgz decompression stream.
  [[nodiscard]] static std::optional<std::unique_ptr<BaseStreamIO>> getStreamIO( const std::string& file_name
                                                                                , size_t decompression_threads = BGZ_DEFAULT_THREADS);


protected:

  // Number of threads used in the decompression pipeline.
  constexpr static size_t BGZ_DEFAULT_THREADS{15};

  constexpr static const char* GZ_FILE_EXTENSION_ = ".GZ"; // gzipped file assumed (checked for '.bgz' format).
  constexpr static const char* BGZ_FILE_EXTENSION_ = ".BGZ"; // gzipped file assumed.
  constexpr static const char* BZ2_FILE_EXTENSION_ = ".BZ2"; // Burrows-Wheeler compression assumed.

};


} // namespace


#endif //KEL_BASIC_IO_H
