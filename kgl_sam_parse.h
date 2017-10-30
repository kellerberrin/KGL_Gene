//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SAM_PARSE_H
#define KGL_SAM_PARSE_H


#include <string>
#include <vector>
#include <queue>
#include "kgl_logging.h"
#include "kgl_mt_data.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class SAMRecordParser {

public:

  explicit SAMRecordParser(Logger& logger) : log(logger) {}
  ~SAMRecordParser() = default;

  // Returns false if the sam record is unmapped.
  bool parseSAM(std::unique_ptr<const std::string>& record_ptr);

  inline const std::vector<std::pair<std::size_t, std::size_t>>& getOptFlags() { return opt_flags_; }

  inline const std::vector<std::pair<const char, const ContigOffset_t>>& cigarFields() { return cigar_fields_; }

  inline const ContigOffset_t getPos(std::unique_ptr<const std::string>& record_ptr) {

    // Subtract 1 for zero offset - SAM usage is offset 1
    return std::stoull(record_ptr->substr(sam_fields_[POS_OFFSET].first, sam_fields_[POS_OFFSET].second)) - 1;

  }
  inline const CharSequence_t getSubSequence( std::unique_ptr<const std::string>& record_ptr
  , std::size_t start
  , std::size_t length) {

    return record_ptr->substr(sam_fields_[SEQUENCE_OFFSET].first + start, length);

  }
  inline const Nucleotide_DNA5_t getSequenceNucleotide(std::unique_ptr<const std::string>& record_ptr, std::size_t offset) {

    return record_ptr->at(sam_fields_[SEQUENCE_OFFSET].first + offset);

  }
  inline const Nucleotide_DNA5_t getQualityNucleotide(std::unique_ptr<const std::string>& record_ptr, std::size_t offset) {

    return record_ptr->at(sam_fields_[QUALITY_OFFSET].first + offset);

  }
  inline const ContigId_t getContigId(std::unique_ptr<const std::string>& record_ptr) {

    return record_ptr->substr(sam_fields_[RNAME_OFFSET].first, sam_fields_[RNAME_OFFSET].second);

  }

private:

  Logger& log;        // Declared First. Emit log messages to console and log file.

  // Define the SAM record offsets.
  static constexpr size_t QNAME_OFFSET = 0;
  static constexpr size_t FLAG_OFFSET = 1;
  static constexpr size_t RNAME_OFFSET = 2;
  static constexpr size_t POS_OFFSET = 3;
  static constexpr size_t MAP_QUALITY_OFFSET = 4;
  static constexpr size_t CIGAR_OFFSET = 5;
  static constexpr size_t RNEXT_OFFSET = 6;
  static constexpr size_t PNEXT_OFFSET = 7;
  static constexpr size_t TLEN_OFFSET = 8;
  static constexpr size_t SEQUENCE_OFFSET = 9;
  static constexpr size_t QUALITY_OFFSET = 10;

  static constexpr int SAM_FIELD_COUNT = 11;
  static constexpr char SAM_DELIMITER = '\t';
  static constexpr const char UNMAPPED_READ_{'*'};

  // Just store the offset and size of fields within the SAM record.
  std::pair<std::size_t, std::size_t> sam_fields_[SAM_FIELD_COUNT];
  std::vector<std::pair<std::size_t, std::size_t>> opt_flags_;
  std::vector<std::pair<const char, const ContigOffset_t>> cigar_fields_;

  bool parseSAMFields(std::unique_ptr<const std::string>& record_ptr, bool parse_opt_fields); // False if unmapped read
  bool fastDecodeSAMCigar( std::unique_ptr<const std::string>& record_ptr // False if unexpected cigar format
  , const std::pair<std::size_t
  , std::size_t>& cigar_offset);

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SAM_PARSE_H
