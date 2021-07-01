//
// Created by kellerberrin on 9/1/21.
//

#include "kgl_variant_db_type.h"

namespace kgl = kellerberrin::genome;


kgl::DataCharacteristic kgl::DataDB::dataCharacteristic() const {

  auto type = findCharacteristic(data_source_);

  if (not type) {

    // Should never happen.
    ExecEnv::log().critical("DataDB::dataCharacteristic; critical error unknown population type, program terminates.");

  }

  return type.value();

}


std::optional<kgl::DataCharacteristic> kgl::DataDB::findCharacteristic(const std::string& source_text)  {

  for (auto const& type : data_characteristics_) {

    if (source_text == type.source_text) return type;

  }

  return std::nullopt;

}

std::optional<kgl::DataCharacteristic> kgl::DataDB::findCharacteristic(DataSourceEnum data_source) {

  for (auto const& type : data_characteristics_) {

    if (data_source == type.data_source) return type;

  }

  return std::nullopt;

}


std::string kgl::DataDB::dataSource(DataSourceEnum data_source) {

  auto source_char_opt = findCharacteristic(data_source);
  if (not source_char_opt) {

    return "";

  }

  return source_char_opt.value().source_text;

}
