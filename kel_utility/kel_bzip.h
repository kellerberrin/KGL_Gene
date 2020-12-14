//
// Created by kellerberrin on 13/12/20.
//

#ifndef KEL_BZIP_H
#define KEL_BZIP_H

#include "kel_mt_queue.h"
#include "kel_basic_io.h"


#include <string>
#include <memory>
#include <future>


namespace kellerberrin {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Multi-threaded decompression of block encoded VCF files "chromosome.vcf.bgz".
// Each block is individually zlib deflated by a thread and the resultant vcf line records
// are queued for further processing.
// This avoids a bottleneck in processing large VCF files as library decompression code only uses
// a single thread to decompress large VCF files.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This struct maps to the gzip header and we need to be careful about any structure byte padding.
// However the layout of the structure and GCC struct padding rules seem to preclude any padding.

struct GZHeaderblock {

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

struct GZTrailerBlock {

  uint32_t crc_check;                     // CRC of decompressed data.
  uint32_t uncompressed_size;             // For compressed VCF files must be in the range [1, 65536]

};

// Returned from the data decompression threads wrapped in a std::future.
// An empty Uncompressed block is used as an eof marker.
struct UncompressedBlock {

  UncompressedBlock() =default;
  ~UncompressedBlock() =default;
  UncompressedBlock(UncompressedBlock&& copy) noexcept {

    block_id = copy.block_id;
    uncompressed_data = std::move(copy.uncompressed_data);
    uncompressed_size = copy.uncompressed_size;

  }

  size_t block_id{0};
  std::unique_ptr<char[]> uncompressed_data;
  size_t uncompressed_size{0};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Decompression object.
// Reads and verifies the compressed .bgz file
// Decompresses the data blocks using a thread pool.
// Enqueues the decompressed (using zlib) data in a thread safe queue.
// Presents a stream-like readLine() interface to the data consumer.
// Data ReadLine() records are guaranteed to be read sequentially and will block
// until the next logical sequential line record is available.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GZBlockDecompression {

public:

  explicit GZBlockDecompression(std::string file_name) : file_name_(std::move(file_name)) {}
  ~GZBlockDecompression() = default;

  // Checks the internal data structures of a VCF bgz file.
  bool verifyGZBlockFile();
  // Decompresses a VCF bgz file using multiple threads and enqueues the decompressed data blocks.
  bool decompressGZBlockFile(size_t thread_count);

private:

  const std::string file_name_;
  // Blocks are queued here to be decompressed.
  // Queue high tide and low tide markers are guessed as reasonable values.
  constexpr static const size_t QUEUE_LOW_TIDE_{100};
  constexpr static const size_t QUEUE_HIGH_TIDE_{500};
  BoundedMtQueue<std::future<UncompressedBlock>> decompress_queue_{QUEUE_HIGH_TIDE_, QUEUE_LOW_TIDE_};
  // Queues parsed line records.
  constexpr static const size_t LINE_LOW_TIDE_{10000};
  constexpr static const size_t LINE_HIGH_TIDE_{50000};
  BoundedMtQueue<IOLineRecord> line_queue_{LINE_HIGH_TIDE_, LINE_LOW_TIDE_};

  // These constants are used to verify the structure of the .bgz file.
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
  constexpr static const size_t HEADER_SIZE_{sizeof(GZHeaderblock)};
  constexpr static const size_t TRAILER_SIZE_{sizeof(GZTrailerBlock)};
  constexpr static const size_t MAX_UNCOMPRESSED_SIZE_{65536};
  constexpr static const size_t BLOCK_SIZE_ADJUST_{HEADER_SIZE_ + TRAILER_SIZE_ - 1};
  constexpr static const size_t INFLATE_WINDOW_FLAG_ = 15 + 32;

  // Thread pool worker function, calls the zlib inflate function.
  [[nodiscard]] static UncompressedBlock decompressBlock( size_t block_count,
                                                          std::byte* compressed_data,
                                                          size_t compressed_data_size);
  // Manage the line queue
  void queueLines();

};



} // Namespace.


#endif //KEL_BZIP_H
