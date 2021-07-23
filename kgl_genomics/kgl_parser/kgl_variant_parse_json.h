//
// Created by kellerberrin on 23/7/21.
//

#ifndef KGL_VARIANT_INFO_JSON_H
#define KGL_VARIANT_INFO_JSON_H


#include "kgl_genome_types.h"
#include "kgl_data_file_impl.h"

#include "kel_bound_queue.h"

#include <memory>
#include <string>
#include <vector>


namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class JSONInfoParser {

public:

  JSONInfoParser() = default;
  ~JSONInfoParser() = default;

  // Begin reading IO records, spawns threads.
  [[nodiscard]] bool commenceJSONIO(const std::string& json_file_name);

  [[nodiscard]] const std::string& getFileName() const { return file_data_.fileName(); }
  [[nodiscard]] size_t dequeueJSONline();

private:

  // The upstream queue of line records.
  FileDataIO file_data_;

  static const constexpr size_t REPORT_INTERVAL{10000}; // Parser progress messages..

};



} // end namespace



#endif // KGL_VARIANT_INFO_JSON_H
