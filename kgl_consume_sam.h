// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
//
//
// Created by kellerberrin on 18/09/17.
//

#ifndef KGL_CONSUME_SAM_H
#define KGL_CONSUME_SAM_H

#include <string>
#include <vector>
#include <queue>
#include "kgl_logging.h"
#include "kgl_mt_data.h"
#include "kgl_parse_sam.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Consume SAM records coming from the SAM record reader (the producer).
// This code is multi-threaded, all data structures must be thread safe (see kgl_mt_data.h).

template <class ConsumerRecordType>
class ConsumeMTSAM {

public:

  explicit ConsumeMTSAM(Logger& logger, std::shared_ptr<ConsumerRecordType> data_map_ptr) : log(logger), contig_data_map_(data_map_ptr)  {}
  virtual ~ConsumeMTSAM() = default;

  ConsumerRecordType& contigDataMap() { return *contig_data_map_; }

  void consume(std::unique_ptr<const std::string>& record_ptr);  // Parse SAM record into fields.
  void finalize();
  void readQuality(unsigned char read_quality) { read_quality_ = read_quality + NUCLEOTIDE_QUALITY_ASCII; }

private:

  Logger& log;                              // Declared First. Emit log messages to console and log file.

  static constexpr unsigned char NUCLEOTIDE_QUALITY_ASCII = 33;  // adjusts nucleotide quality to an ascii value.
  static constexpr unsigned char DEFAULT_MINIMUM_QUALITY = 30 + NUCLEOTIDE_QUALITY_ASCII;   // -10 log10 Pr{ReadError}

  unsigned char read_quality_ = DEFAULT_MINIMUM_QUALITY;
  std::shared_ptr<ConsumerRecordType> contig_data_map_;             // Thread safe map of contig data.

  std::atomic<uint64_t> unmapped_reads_{0};     // Contig = "*"
  std::atomic<uint64_t> other_contig_{0};       // SAM records for unregistered contigs (ignored by the code).
  std::atomic<uint64_t> insert_sequence_{0};    // Insertions
  std::atomic<uint64_t> delete_nucleotide_{0};  // Deletions
  std::atomic<uint64_t> rejected_{0};           // Below Quality Reads
  std::atomic<uint64_t> accepted_{0};           // Nucleotides accepted

};


template <class ConsumerRecordType>
void ConsumeMTSAM<ConsumerRecordType>::finalize() {

  size_t unmapped = unmapped_reads_;
  size_t deleted = delete_nucleotide_;
  size_t inserted = insert_sequence_;
  size_t ignored = other_contig_;

  log.info("Processed SAM file; unmapped: {}, deleted: {}, inserted: {}, ignored: {}"
      , unmapped, deleted, inserted, ignored);

  if (read_quality_ > NUCLEOTIDE_QUALITY_ASCII) {

    size_t rejected = rejected_;
    size_t accepted = accepted_;
    double percentage = static_cast<double>(rejected)/static_cast<double>(accepted) * 100.0;

    log.info("Accept/Reject nucleotide read quality -10log10(Pr Read Error) = {}, accepted {}, rejected: {} ({}%)",
             static_cast<int>(read_quality_ - NUCLEOTIDE_QUALITY_ASCII), accepted , rejected, percentage);

  }
  else {

    log.info("Accept/Reject nucleotide read quality disabled (read quality = 0)");

  }

}

template <class ConsumerRecordType>
void ConsumeMTSAM<ConsumerRecordType>::consume(std::unique_ptr<const std::string>& record_ptr) {

  SAMRecordParser sam_record_parser(log);

  if (not sam_record_parser.parseSAM(record_ptr)) {

    ++unmapped_reads_;  // Contig is "*"
    return;

  }

  const ContigId_t contig_id(sam_record_parser.getContigId(record_ptr));

  auto contig_data_ptr = contig_data_map_->findContigBlock(contig_id);

  if (not contig_data_ptr) {

    ++other_contig_;    // Contig is not registered.
    return;

  }

  auto& contig_block = contig_data_ptr->getNucleotideArray();
  auto contig_size = contig_block.contigSize();
  auto& insert_block = contig_data_ptr->getInsertArray();

  ContigOffset_t location = sam_record_parser.getPos(record_ptr);

  if (location >= contig_size) {

    log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {}; SAM record: {}",
              contig_id, contig_size, location, *record_ptr);
    return;

  }

  ContigOffset_t sam_idx = 0;

  for (auto cigar : sam_record_parser.cigarFields()) {

    switch(cigar.first) {

      case 'M':
      case 'X':
      case '=':
        for (ContigOffset_t cigar_offset = 0; cigar_offset < cigar.second; ++cigar_offset) {

          if (location + cigar_offset >= contig_size) {
            log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {} Cigar Offset: {}; SAM record: {}"
                , contig_id, contig_size, location, cigar_offset, *record_ptr);
            return;
          }
          if (read_quality_ > NUCLEOTIDE_QUALITY_ASCII) {

            if (sam_record_parser.getQualityNucleotide(record_ptr, cigar_offset + sam_idx) >= read_quality_) {

              Nucleotide_t sam_nucleotide = sam_record_parser.getSequenceNucleotide(record_ptr, cigar_offset + sam_idx);
              contig_block.incrementCount(location + cigar_offset, sam_nucleotide);
              ++accepted_;

            } else {

              ++rejected_;

            }

          } else { // Read quality disabled.

            Nucleotide_t sam_nucleotide = sam_record_parser.getSequenceNucleotide(record_ptr, cigar_offset + sam_idx);
            contig_block.incrementCount(location + cigar_offset, sam_nucleotide);

          }

        }

        sam_idx += cigar.second;
        location += cigar.second;
        break;

      case 'D':

        ++delete_nucleotide_;
        for (ContigOffset_t cigar_offset = 0; cigar_offset < cigar.second; ++cigar_offset) {

          contig_block.incrementDelete(location + cigar_offset);

        }

        location += cigar.second;
        break;

      case 'I':

        ++insert_sequence_;
        contig_block.incrementInsert(location);
        insert_block.insertSequence(location, sam_record_parser.getSubSequence(record_ptr, sam_idx, cigar.second));

        sam_idx += cigar.second;
        break;

      case 'S':
        sam_idx += cigar.second;
        break;

      case 'H':
      case 'P':
        break;

      case 'N':
        location += cigar.second;
        break;

      default:
        log.error("Unknown cigar code {}; SAM record: {}", cigar.first, *record_ptr);
        return;

    }

  }

}


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_CONSUME_SAM_H
