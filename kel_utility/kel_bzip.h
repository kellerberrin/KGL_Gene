//
// Created by kellerberrin on 13/12/20.
//

#ifndef KEL_BZIP_H
#define KEL_BZIP_H

#include "kel_queue_tidal.h"
#include "kel_workflow_threads.h"
#include "kel_basic_io.h"

#include <string>
#include <memory>


namespace kellerberrin {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// With 15 decompression threads, Will decompress a 200gb '.bgz' file to 1.45tb
// of individual line records (63,000,000) in about 300 seconds.
// The use of 'tidal' queues to hold intermediate results means that memory usage is insignificant (< 10mb).
//
// The block gzip file format decompressed by this code is specified by standard RFC1952.
// This is the format used in compressing large genome VCF files.
// This code will not decompress standard .gz files.
// The Block gzip (.bgz) decompression object below reads and verifies compressed .bgz (RFC1952) files.
// The file is decompressed the using a thread-pool.
// This avoids a bottleneck in processing large VCF files as library decompression code only uses
// a single thread to decompress large VCF files.
// The workflow enqueues the decompressed (using zlib) data in a thread safe queue.
// The object presents a stream-like readLine() interface to the data consumer.
// Data ReadLine() records are guaranteed to be read sequentially and will block
// until the next logical sequential line record is available (except on eof which does not block).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Current state of the BGZReader, 'ACTIVE' and processing after open() or 'STOPPED' before open() or after close().
// Note that the BGZReader is still 'ACTIVE' on an eof condition. It is only 'STOPPED' after calling close().
// The BGZReader can read multiple '.bgz' files by sequentially calling close() and then open().
enum class BGZReaderState { ACTIVE, STOPPED};

class BGZReader : public BaseStreamIO {

  // This struct maps to the byte structure of the '.bgz' header and we need to be careful about any byte padding.
  // The layout of the structure and GCC struct padding rules seem to preclude any padding and so this appears to be OK.
  // But this could be a source of grief with another compiler with different byte padding rules to GCC.
  // Hopefully the alignas(alignof(...)) specifiers below will suppress any compiler specific structure padding.
  struct BGZHeaderblock {

    alignas(alignof(uint8_t)) uint8_t block_id_1;                // aligned on a single byte, value 31
    alignas(alignof(uint8_t)) uint8_t block_id_2;                // aligned on a single byte, value 139
    alignas(alignof(uint8_t)) uint8_t compression_method;        // aligned on a single byte, value 8
    alignas(alignof(uint8_t)) uint8_t flags;                     // aligned on a single byte, value 4.
    alignas(alignof(uint32_t)) uint32_t mtime;                   // aligned on 4 byte boundary, OK because preceded by 4 single bytes.
    alignas(alignof(uint8_t)) uint8_t extra_flags;               // aligned on a single byte
    alignas(alignof(uint8_t)) uint8_t operating_system;          // aligned on a single byte
    alignas(alignof(uint16_t)) uint16_t length_extra_blocks;     // aligned on 2 byte boundary, OK, value is 6.
    alignas(alignof(uint8_t)) uint8_t subfield_id_1;             // aligned on a single byte, value 66
    alignas(alignof(uint8_t)) uint8_t subfield_id_2;             // aligned on a single byte, value 67
    alignas(alignof(uint16_t)) uint16_t subfield_length;         // aligned on 2 byte boundary, OK, value 2.
    // Total Block SIZE (include header + trailer) minus 1.
    alignas(alignof(uint16_t)) uint16_t block_size;              // aligned on 2 byte boundary,
    // Variable compressed data size is block_size - (header_size + trailer_size - 1)

  };

  struct BGZTrailerBlock {

    alignas(alignof(uint32_t)) uint32_t crc_check;                // CRC of decompressed data.
    alignas(alignof(uint32_t)) uint32_t uncompressed_size;        // For compressed VCF files must be in the range [1, 65536]

  };

  // Holds the compressed data read directly from file.
  // Note the union to hold the header block.
  constexpr static const size_t MAX_UNCOMPRESSED_SIZE_{65536};
  using BGZCompressedData = std::array<std::byte, MAX_UNCOMPRESSED_SIZE_>;
  struct CompressedBlock {

    union {
      alignas(alignof(uint32_t)) BGZHeaderblock header_block_;
      alignas(alignof(uint32_t)) BGZCompressedData compressed_block_;
    };
    size_t block_id_{0};
    size_t data_size_{0};

  };

  // Returned from the data decompression threads wrapped in a std::future.
  struct DecompressedBlock {

    std::vector<std::unique_ptr<std::string>> parsed_records;
    size_t block_id{0};

  };

  // Convenience typedefs.
  using CompressedType = std::unique_ptr<CompressedBlock>;
  using DecompressedType = std::unique_ptr<DecompressedBlock>;

public:

  explicit BGZReader(size_t decompression_threads = DECOMPRESSION_THREADS_) : decompress_threads_(decompression_threads) {}
  ~BGZReader() override { close(); }

  // Guaranteed sequential line reader. Does not block on EOF.
  [[nodiscard]] IOLineRecord readLine() override;

  // Opens the '.bgz' file and begins decompressing the file, on success object is now 'ACTIVE'.
  bool open(const std::string &file_name) override;

  // Close the physical file and reset the internal queues and threads, the object is now 'STOPPED'.
  // Once the object is closed, it may be re-opened to process another file.
  bool close();

  [[nodiscard]] bool good() const { return bgz_file_.good(); }

  // Checks the internal data structures of a RFC1952 standard '.bgz' file. May be slow on very large files.
  // Returns true for a valid '.bgz' with a RFC1952 structure. A normal (non RFC1952) '.gz' will return false.
  [[nodiscard]] static bool verify(const std::string &file_name, bool silent = true);

  // State of the object, 'ACTIVE' or 'STOPPED'.
  [[nodiscard]] BGZReaderState readerState() const { return reader_state_; }

  // Probes to examine the queues for optimal tidal size and thread allocation.
  [[nodiscard]] const QueueTidal<std::future<DecompressedType>>& decompressQueue() const { return decompress_queue_; }
  [[nodiscard]] const QueueTidal<IOLineRecord>& lineQueue() const { return line_queue_; }

private:

  // Number of threads used in the decompression pipeline.
  constexpr static const size_t DECOMPRESSION_THREADS_{15};

  std::string file_name_;
  std::ifstream bgz_file_;
  WorkflowThreads decompress_threads_;   // Multiple threads to decompress the bgz data blocks.
  WorkflowThreads reader_thread_{1}; // Reads from the file and passes compressed bgz blocks to the decompression threads.
  WorkflowThreads line_asssemble_thread_{1};  // Assembles complete line records from the decompressed data.
  std::atomic<bool> close_stream_{false};  // Set to gracefully shutdown the decompression pipeline when close() is called.
  std::atomic<BGZReaderState> reader_state_{BGZReaderState::STOPPED};
  bool line_eof_{false};  // Set on the file EOF condition when there are no more blocks to be decompressed.

  // 'Tidal' (bounded) queue for compressed blocks to be decompressed.
  // Queue high tide and low tide markers are tested as reasonable values.
  constexpr static const size_t QUEUE_LOW_TIDE_{2000};
  constexpr static const size_t QUEUE_HIGH_TIDE_{4000};
  constexpr static const char* QUEUE_NAME_{"BGZReader Decompress Block Queue"};
  constexpr static const size_t QUEUE_SAMPLE_FREQ_{500};
  QueueTidal<std::future<DecompressedType>> decompress_queue_{QUEUE_HIGH_TIDE_, QUEUE_LOW_TIDE_, QUEUE_NAME_, QUEUE_SAMPLE_FREQ_};

  // 'Tidal' (bounded) queue for parsed line records.
  constexpr static const size_t LINE_LOW_TIDE_{10000};
  constexpr static const size_t LINE_HIGH_TIDE_{20000};
  constexpr static const char* LINE_QUEUE_NAME_{"BGZReader Line Record Queue"};
  constexpr static const size_t LINE_SAMPLE_FREQ_{500};
  QueueTidal<IOLineRecord> line_queue_{LINE_HIGH_TIDE_, LINE_LOW_TIDE_, LINE_QUEUE_NAME_, LINE_SAMPLE_FREQ_};

  // These constants are used to verify the structure of the '.bgz' file.
  // Don't change these constants, they are specified by standard RFC1952.
  constexpr static const uint8_t BLOCK_ID1_{31};
  constexpr static const uint8_t BLOCK_ID2_{139};
  constexpr static const uint8_t COMPRESSION_{8};  // 8 = Z_DEFLATED in zlib.h
  constexpr static const uint8_t FLAGS_{4};
  constexpr static const uint8_t EXTRA_LENGTH_{6};
  constexpr static const uint8_t SUBFIELD_ID1_{66};
  constexpr static const uint8_t SUBFIELD_ID2_{67};
  constexpr static const uint8_t SUBFIELD_LENGTH_{2};

  constexpr static const size_t EOF_MARKER_SIZE_{28};
  constexpr static const uint8_t EOF_MARKER_[28] { 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
                                                   0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
                                                   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
  constexpr static const size_t HEADER_SIZE_{sizeof(BGZHeaderblock)};
  constexpr static const size_t TRAILER_SIZE_{sizeof(BGZTrailerBlock)};
  constexpr static const size_t BLOCK_SIZE_ADJUST_{HEADER_SIZE_ + TRAILER_SIZE_ - 1};
  constexpr static const size_t INFLATE_WINDOW_FLAG_ = 15 + 32;
  constexpr static const char EOL_MARKER_ = '\n';

  // Reads and enqueues the compressed data blocks.
  void decompressGZBlockFile();
  // Thread pool worker function, calls the zlib inflate function. Most work is done here.
  [[nodiscard]] static DecompressedType decompressBlock(std::shared_ptr<CompressedBlock> compressed_data);
  // Assemble records last + first and queue as complete records.
  void assembleRecords();
  // Check the trailing '.bgz' file EOF_MARKER_.
  void checkEOFMarker(size_t remaining_chars);

};



} // Namespace.


#endif //KEL_BZIP_H
