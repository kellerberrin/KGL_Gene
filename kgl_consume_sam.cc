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
// Created by kellerberrin on 18/09/17.
//

#include "kgl_parse_sam.h"
#include "kgl_consume_sam.h"

namespace kgl = kellerberrin::genome;


void kgl::ConsumeMTSAM::finalize() {

  size_t unmapped = unmapped_reads_;
  size_t deleted = delete_nucleotide_;
  size_t inserted = insert_sequence_;
  size_t ignored = other_contig_;
  log.info("Completed processing SAM file; unmapped: {}, deleted: {}, inserted: {}, ignored: {}"
      , unmapped, deleted, inserted, ignored);

}

void kgl::ConsumeMTSAM::consume(std::unique_ptr<const std::string>& record_ptr) {

  SAMRecordParser sam_record_parser(log);

  if (not sam_record_parser.parseSAM(record_ptr)) {

    ++unmapped_reads_;  // Contig is "*"
    return;

  }

  const ContigId_t contig_id(sam_record_parser.getContigId(record_ptr));

  auto matrix_ptr = contig_data_map_.getContig(contig_id);  // Get the contig data block.

  if (matrix_ptr == contig_data_map_.notFound()) {

    ++other_contig_;    // Contig is not registered.
    return;

  }

  auto& contig_block = contig_data_map_.getMatrix(matrix_ptr).getNucleotideArray();
  auto& insert_block = contig_data_map_.getMatrix(matrix_ptr).getInsertArray();

  auto contig_size = contig_block.contigSize();

  ContigOffset_t location = sam_record_parser.getPos(record_ptr);

  if (location >= contig_size) {

    log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {}; SAM record: {}"
        , contig_id, contig_size, location, *record_ptr);
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
          Nucleotide_t sam_nucleotide = sam_record_parser.getSequenceNucleotide(record_ptr, cigar_offset + sam_idx);
          contig_block.incrementCount(location + cigar_offset, sam_nucleotide);

        }

        sam_idx += cigar.second;
        location += cigar.second;
        break;

      case 'D':

        delete_nucleotide_++;
        for (ContigOffset_t cigar_offset = 0; cigar_offset < cigar.second; ++cigar_offset) {

          contig_block.incrementDelete(location + cigar_offset);

        }

        location += cigar.second;
        break;

      case 'I':

        insert_sequence_++;
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

