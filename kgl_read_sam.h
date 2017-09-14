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
// Created by kellerberrin on 8/09/17.
//

#ifndef KGL_READ_SAM_H
#define KGL_READ_SAM_H


#include <string>
#include <vector>
#include <queue>
#include "kgl_logging.h"
#include "kgl_mt_queue.h"
#include "kgl_process_sam.h"
#include "kgl_mt_data.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ProcessSamFile {

public:

  explicit ProcessSamFile(const std::string& log_file);
  virtual ~ProcessSamFile() = default;  // Called by the Python binding super class.

  void readSamFile(std::string &file_name);   // Mainline spawns the consumer threads.
  ContigDataMap& contigDataMap() { return contig_data_map_; }
  InsertQueue& getInsertQueue() { return  insert_queue_ ; };


private:

  static constexpr const char* sam_read_module_name_{"SamRead"};  // Name of this module for the logger
  static constexpr const char* eof_indicator_{"<<EOF>>"};  // Enqueued by producer to indicate SAM eof.
  static constexpr long report_increment_{500000};    // Frequency to emit SAM progress messages
  static constexpr long high_tide_{1000000};          // Maximum BoundedMtQueue size
  static constexpr long low_tide_{500000};            // Low water mark to begin queueing SAM records
  static constexpr char delete_nucleotide{'-'};
  static constexpr char insert_nucleotide{'+'};
  static constexpr const char* unmapped_read{"*"};

  Logger log;                                         // Emit log messages to console and log file.

  int consumer_thread_count_{4};                      // Consumer threads (defaults to local CPU cores available)
  BoundedMtQueue<std::unique_ptr<const std::string>> producer_consumer_queue_; // The Producer/Consumer SAM record queue
  ProcessSamRecord process_sam_record_;               // Must be declared after the logger.
  ContigDataMap contig_data_map_;                     // Must be declared after the logger.
  InsertQueue insert_queue_;       // Enqueued by spawned threads and dequeued by mainline.

  std::atomic<uint64_t> unmapped_reads_{0};  // Read statistics.
  std::atomic<uint64_t> insert_sequence_{0};
  std::atomic<uint64_t> delete_nucleotide_{0};
  std::atomic<uint64_t> mismatch_nucleotide_{0};

  void samProducer(std::string &file_name);   // Read the SAM file and queue the record in a BoundedMtQueue.
  void samConsumer();                         // Multiple threads; dequeue from the BoundedMtQueue and process.
  void parseSAMRecord(std::unique_ptr<const std::string>& record_ptr);  // Parse SAM record into fields.
  void parseSAMFields( std::string sam_record
                     , std::vector<std::string>& sam_fields
                     , std::vector<std::string>& sam_flags);   // Split record into fields
  void decodeSAMCigar( std::string cigar_string
                     , std::vector<std::pair<const char, const ContigOffset_t >>& cigar_fields);

};



}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_READ_SAM_H
