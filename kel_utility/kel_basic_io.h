//
// Created by kellerberrin on 11/5/20.
//

#ifndef KEL_BASIC_IO_H
#define KEL_BASIC_IO_H

#include "kel_queue_tidal.h"
#include "kel_workflow_threads.h"

#include <string>
#include <string_view>
#include <vector>
#include <optional>
#include <fstream>


namespace kellerberrin {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A lightweight object that contains the readLine() output of StreamIO with file line number and associated line text.
// Note that the line text has been stripped of the terminating '\n'.
// Once the line data has been moved from the object using getLineData() it is tagged as empty().
// A subsequent attempt to access the moved data using getLineData() will generate a warning.
// A static function createEOFMarker() creates an object to serve as EOF, this can also be pushed onto a queue.
// Since the stored std::string data can be large, for optimal performance the object cannot be copied.
// The object is designed to be pushed and poped from thread safe queues without incurring significant processing overhead.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class IOLineRecord {

public:

  IOLineRecord(size_t line_count, std::string&& line_data) noexcept : line_count_(line_count), line_data_(std::move(line_data)) {}
  IOLineRecord(IOLineRecord&& line_record) noexcept : line_count_(line_record.line_count_)
                                                    , line_data_(std::move(line_record.line_data_))
                                                    , EOF_(line_record.EOF_)
                                                    , empty_(line_record.empty_) {}
  IOLineRecord(const IOLineRecord& line_record) = delete;
  ~IOLineRecord() = default;

  IOLineRecord& operator=(const IOLineRecord&) = delete;

  // After moving data out of the object, it is empty() but not EOF().
  [[nodiscard]] std::pair<size_t, std::string> getLineData();
  [[nodiscard]] size_t lineCount() const { return line_count_; }
  [[nodiscard]] const std::string_view getView() const { return {line_data_}; }
  [[nodiscard]] bool EOFRecord() const { return EOF_; }
  [[nodiscard]] bool empty() const { return empty_; }


  // Return the EOF marker.
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

class BaseStreamIO {

public:

  BaseStreamIO() = default;
  virtual ~BaseStreamIO() = default;

  virtual bool open(const std::string &file_name) = 0;
  virtual IOLineRecord readLine() = 0; // This function is NOT multithreaded.

  // Returns an IO stream that is either a normal stream or a compressed stream based on the file name extension.
  // If '.bgz' then it uses a multithreaded and memory efficient algorithm to decompress as it reads.
  // Very large '.bgz' files can be decompressed and streamed to the relevant parser for further analysis.
  // If '.gz' then the file is analyzed to see if it is in '.bgz' format.
  // If not in '.bgz' format then a standard (boost) single thread decompression algorithm is used.
  // If '.bz2' then Burrows-Wheeler compression is assumed.
  // If the file is not one of the above types, it is assumed to be an uncompressed record based text file.
  // Note the stream is returned open and ready for processing, std::nullopt is returned if there is a problem.
  // The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
  [[nodiscard]] static std::optional<std::unique_ptr<BaseStreamIO>> getStreamIO( const std::string& file_name
                                                                               , size_t decompression_threads = BGZ_DEFAULT_THREADS);

  // Default number of threads used in the '.bgz' decompression pipeline.
  constexpr static const size_t BGZ_DEFAULT_THREADS{15};

protected:

  size_t record_counter_{0};

  constexpr static const char* GZ_FILE_EXTENSTION_ = ".GZ"; // gzipped file assumed (checked for '.bgz' format).
  constexpr static const char* BGZ_FILE_EXTENSTION_ = ".BGZ"; // gzipped file assumed.
  constexpr static const char* BZ2_FILE_EXTENSTION_ = ".BZ2"; // Burrows-Wheeler compression assumed.

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A queued multi-thread adapter for text Stream IO.
// A soon as the underlying stream is successfully opened, records are read and stored in a tidal queue and can be
// retrieved using readLine().
// This object cannot be copied.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MTStreamIO {

public:

  MTStreamIO() = default;
  ~MTStreamIO() = default;

  // Uses the filename extension heuristic documented above to open the underlying StreamIO
  // The threads argument is only valid for '.bgz' file types. The argument is ignored for other stream types.
  bool open(const std::string &file_name, size_t decompression_threads = BaseStreamIO::BGZ_DEFAULT_THREADS);
  // After closing the object can be re-opened on another file.
  void close();

  // This is thread safe.
  // If read by 1 thread then records are guaranteed to presented in file order.
  // If read by multiple threads then there is no line order guarantee.
  // This function will block on an empty queue and no EOF.
  // After an EOF has been received, then subsequent calls will not block and will return EOF objects.
  [[nodiscard]] IOLineRecord readLine();
  // Access queue stats.
  [[nodiscard]] const QueueTidal<IOLineRecord>& lineQueue() const { return line_io_queue_; }

private:

  // The tidal IO queue parameters.
  static constexpr const size_t IO_HIGH_TIDE_{10000};          // Maximum QueueTidal size
  static constexpr const size_t IO_LOW_TIDE_{2000};            // Low water mark to begin queueing data records
  static constexpr const char* IO_QUEUE_NAME_{"MTStreamIO Queue"};      // The queue name
  static constexpr const size_t IO_SAMPLE_RATE_{100};            // The queue monitor sampling rate (ms), zero (0) disables the monitor.
  QueueTidal<IOLineRecord> line_io_queue_{IO_HIGH_TIDE_, IO_LOW_TIDE_, IO_QUEUE_NAME_, IO_SAMPLE_RATE_};
  // Queues IOLineRecord objects to the queue.
  static constexpr const size_t WORKER_THREAD_COUNT{1};
  WorkflowThreads line_io_thread_;
  // The StreamIO object.
  std::unique_ptr<BaseStreamIO> stream_ptr_;
  // Set if the stream is in an EOF condition.
  std::atomic<bool> EOF_received_{false};
  // Mutex for readLine()
  std::mutex mutex_;
  // The thread worker function.
  void enqueueIOLineRecord();

};


} // namespace


#endif //KEL_BASIC_IO_H
