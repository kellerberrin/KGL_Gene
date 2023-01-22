//
// Created by kellerberrin on 13/12/20.
//

#ifndef KEL_BZIP_H
#define KEL_BZIP_H

#include "kel_bound_queue.h"
#include "kel_thread_pool.h"
#include "kel_basic_io.h"


#include <string>
#include <memory>


namespace kellerberrin {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
// until the next logical sequential line record is available (except on eof).
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class BGZReader : public BaseStreamIO {

  // This struct maps to the byte structure of the gzip header and we need to be careful about any structure byte padding.
  // The layout of the structure and GCC struct padding rules seem to preclude any padding and so this appears to be OK.
  // But this could be a source of grief with another compiler with different byte padding rules to GCC.

  struct BGZHeaderblock {

    uint8_t block_id_1;                       // aligned on a single byte, value 31
    uint8_t block_id_2;                       // aligned on a single byte, value 139
    uint8_t compression_method;               // aligned on a single byte, value 8
    uint8_t flags;                            // aligned on a single byte, value 4.
    uint32_t mtime;                           // aligned on 4 byte boundary, OK because preceded by 4 single bytes.
    uint8_t extra_flags;                      // aligned on a single byte
    uint8_t operating_system;                 // aligned on a single byte
    uint16_t length_extra_blocks;             // aligned on 2 byte boundary, OK, value is 6.
    uint8_t subfield_id_1;                    // aligned on a single byte, value 66
    uint8_t subfield_id_2;                    // aligned on a single byte, value 67
    uint16_t subfield_length;                 // aligned on 2 byte boundary, OK, value 2.
    // Total Block SIZE (include header + trailer) minus 1.
    uint16_t block_size;                      // aligned on 2 byte boundary,
    // Variable compressed data size is block_size - (header_size + trailer_size - 1)

  };

  struct BGZTrailerBlock {

    uint32_t crc_check;                     // CRC of decompressed data.
    uint32_t uncompressed_size;             // For compressed VCF files must be in the range [1, 65536]

  };

  // Holds the compressed read directly from file. Note the union to hold the header block.
  constexpr static const size_t MAX_UNCOMPRESSED_SIZE_{65536};
  using BGZCompressedData = std::array<std::byte, MAX_UNCOMPRESSED_SIZE_>;
  struct CompressedBlock {

    union {
      BGZHeaderblock header_block_;
      BGZCompressedData compressed_block_;
    };
    size_t block_id_{0};
    size_t data_size_{0};

  };

  // Returned from the data decompression threads wrapped in a std::future.
  struct UncompressedBlock {

    std::vector<std::unique_ptr<std::string>> parsed_records;
    size_t block_id{0};

  };

  // Convenience alias.
  using CompressedType = std::unique_ptr<CompressedBlock>;
  using DecompressedType = std::unique_ptr<UncompressedBlock>;

public:

  explicit BGZReader(size_t thread_count = DEFAULT_THREADS) : decompress_threads_(thread_count) {}
  ~BGZReader() override { close(); }

  // Guaranteed sequential line reader. Does not block on EOF.
  IOLineRecord readLine() override;

  bool open(const std::string &file_name) override;

  bool close();

  bool good() const { return bgz_file_.good(); }

  // Checks the internal data structures of a VCF bgz file.
  static bool verify(const std::string &file_name, bool silent = true);

  // Seems about right.
  constexpr static const size_t DEFAULT_THREADS{15};

  // Examine the queues for optimal tidal size and thread allocation.
  const BoundedMtQueue<std::future<DecompressedType>>& decompressQueue() const { return decompress_queue_; }
  const BoundedMtQueue<IOLineRecord>& lineQueue() const { return line_queue_; }

private:

  std::string file_name_;
  std::ifstream bgz_file_;
  WorkflowThreads decompress_threads_;   // Multiple threads to decompress the bgz data blocks.
  WorkflowThreads reader_thread_{1}; // Reads from the file and passes compressed bgz blocks to the decompression threads.
  WorkflowThreads line_asssemble_thread_{1};  // Assembles complete line records from the decompressed data.
  bool line_eof_{false};

  // Blocks are queued here to be decompressed.
  // Queue high tide and low tide markers are tested as reasonable values.
  constexpr static const size_t QUEUE_LOW_TIDE_{2000};
  constexpr static const size_t QUEUE_HIGH_TIDE_{4000};
  constexpr static const char* QUEUE_NAME_{"BGZReader Decompress Block Queue"};
  constexpr static const size_t QUEUE_SAMPLE_FREQ_{500};
  BoundedMtQueue<std::future<DecompressedType>> decompress_queue_{QUEUE_HIGH_TIDE_, QUEUE_LOW_TIDE_, QUEUE_NAME_, QUEUE_SAMPLE_FREQ_};

  // Queues parsed line records.
  constexpr static const size_t LINE_LOW_TIDE_{10000};
  constexpr static const size_t LINE_HIGH_TIDE_{20000};
  constexpr static const char* LINE_QUEUE_NAME_{"BGZReader Line Record Queue"};
  constexpr static const size_t LINE_SAMPLE_FREQ_{500};
  BoundedMtQueue<IOLineRecord> line_queue_{LINE_HIGH_TIDE_, LINE_LOW_TIDE_, LINE_QUEUE_NAME_, LINE_SAMPLE_FREQ_};

  // These constants are used to verify the structure of the .bgz file.
  // Don't change these constants
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
  // Thread pool worker function, calls the zlib inflate function.
  [[nodiscard]] static DecompressedType decompressBlock(std::shared_ptr<CompressedBlock> compressed_data);
  // Assemble records last + first and queue as complete records.
  void assembleRecords();
  // Check the trailing EOF_MARKER_
  void checkEOFMarker(size_t remaining_chars);

};



} // Namespace.


#endif //KEL_BZIP_H
