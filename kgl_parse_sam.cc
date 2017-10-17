///
// Created by kellerberrin on 20/09/17.
//

#include <regex>
#include "kgl_parse_sam.h"

namespace kgl = kellerberrin::genome;

// Return false if an unmapped read or problems with cigar.
bool kgl::SAMRecordParser::parseSAM(std::unique_ptr<const std::string>& record_ptr) {

  if (parseSAMFields(record_ptr, false)) {

    if (fastDecodeSAMCigar(record_ptr, sam_fields_[CIGAR_OFFSET])) {

      return true;

    }

  }

  return false;

}

// Efficiently store field offsets and sizes.
bool kgl::SAMRecordParser::parseSAMFields(std::unique_ptr<const std::string>& record_ptr, bool parse_opt_fields)
{

  std::size_t start = 0;
  std::size_t end = record_ptr->find_first_of(SAM_DELIMITER);

  int field_count = 0;

  while (end <= std::string::npos) {

    if (field_count < SAM_FIELD_COUNT) {

      sam_fields_[field_count] = std::pair<std::size_t, std::size_t>(start, end - start);

    } else if (parse_opt_fields) {

      opt_flags_.emplace_back(start, end - start);

    }
    else break;

    ++field_count;

    if (end == std::string::npos)
      break;

    start = end + 1;
    end = record_ptr->find_first_of(SAM_DELIMITER, start);

  }

  if (field_count <  SAM_FIELD_COUNT) {

    log.error("Incorrect field count: {}, SAM record: {}", field_count, *record_ptr);
    return false;

  }

  return (record_ptr->at(sam_fields_[RNAME_OFFSET].first) != UNMAPPED_READ_);

}


bool kgl::SAMRecordParser::fastDecodeSAMCigar( std::unique_ptr<const std::string>& record_ptr
    , const std::pair<std::size_t, std::size_t>& cigar_offset) {

  char cigar_buffer[1000]; // Should big enough
  std::size_t char_index = 0;
  std::size_t buffer_index = 0;

  while (char_index < cigar_offset.second) {

    if (!isdigit(record_ptr->at(cigar_offset.first + char_index))) {

      if (record_ptr->at(cigar_offset.first + char_index) != UNMAPPED_READ_) {

        log.error("Unexpected character in cigar, SAM record: {}", *record_ptr);

      }

      return false;

    }

    while (true) {

      cigar_buffer[buffer_index] = record_ptr->at(cigar_offset.first + char_index);
      if (isalpha(cigar_buffer[buffer_index])) break;
      ++char_index;
      ++buffer_index;
      if (char_index >= cigar_offset.second) {

        return false;

      }

    }
    const char cigar_code = cigar_buffer[buffer_index];
    switch(cigar_code) {

      case 'M':
      case 'X':
      case '=':
      case 'D':
      case 'I':
      case 'S':
      case 'H':
      case 'P':
      case 'N':
        break;

      case UNMAPPED_READ_:
        return false;

      default:
        log.error("Unexpected character in cigar, SAM record: {}", *record_ptr);
        return false;

    }
    cigar_buffer[buffer_index] = '\0';
    ContigOffset_t cigar_length = std::strtoull(cigar_buffer, nullptr, 10);
    cigar_fields_.emplace_back(cigar_code, cigar_length);
    buffer_index = 0;
    ++char_index;

  }

  return true;

}
