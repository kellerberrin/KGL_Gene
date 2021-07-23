//
// Created by kellerberrin on 23/7/21.
//

#include "kgl_variant_parse_json.h"


namespace kgl = kellerberrin::genome;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::JSONInfoParser::commenceJSONIO(const std::string& json_file_name) {

  if (not file_data_.commenceIO(json_file_name)) {

    ExecEnv::log().error("JSONInfoParser::commenceJSONIO; problem processing JSON file: {}", json_file_name);
    return false;

  }

  return true;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse the line record into a VCF record and enqueue the record.

size_t kgl::JSONInfoParser::dequeueJSONline() {

  size_t chars_processed{0};
  size_t records_processed{0};

  while (true) {

    IOLineRecord line_record = file_data_.readIORecord();
    if (not line_record) { // check for EOF condition.

      // push the eof marker back on the queue.
      file_data_.enqueueEOF();
      break;

    }

    auto const& [line_count, line_string_ptr] = line_record.value();

    chars_processed += line_string_ptr->size();
    ++records_processed;

    if (records_processed % REPORT_INTERVAL == 0) {

      ExecEnv::log().info("JSONInfoParser processed line records: {}, line count: {}, text bytes: {}",
                          records_processed, line_count, chars_processed);

    }

  }

  return chars_processed;

}
