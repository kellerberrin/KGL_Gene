//
// Created by kellerberrin on 2/8/20.
//

#include "kgl_ped_parser.h"
#include "kel_utility.h"

namespace kgl = kellerberrin::genome;


void kgl::ParsePedFile::readParsePEDImpl(const std::string& file_name) {


  ExecEnv::log().info("Start Parsing PED file: {}", file_name);

  commenceIO(file_name);

  while (true) {

    IOLineRecord line_record = readIORecord();

    if (not line_record) { // check for EOF condition.

      break;

    }

    auto line_number = line_record.value().first;
    // skip number
    if (line_number == 1) {

      continue;

    }


    if (not moveToPEDRecord(std::move(*line_record.value().second))) {

      ExecEnv::log().warn("ParsePedFile::readParsePEDImpl, Failed to parse PED file: {} record line : {}", fileName(), line_record.value().first);

    }

  }

}


bool kgl::ParsePedFile::moveToPEDRecord(std::string&& line_record) {

  std::vector<std::string> field_strings = Utility::char_tokenizer(line_record, PED_FIELD_DELIMITER_CHAR_);

  if (field_strings.size() != PEDRecord::PEDFieldCount()) {

    ExecEnv::log().error("ParsePedFile::moveToPEDRecord; PED file: {}, record: {} field count: {} not equal mandatory count: {}",
                         fileName(), line_record, field_strings.size(), PEDRecord::PEDFieldCount());
    return false;

  }

  PEDRecord ped_record(field_strings[0],   // family_id
                       field_strings[1],  //individual_id
                       field_strings[2],  // paternal_id
                       field_strings[3],  // maternal_id
                       field_strings[4],  // sex
                       field_strings[5],  // pheno_type
                       field_strings[6],  // population
                       field_strings[7],  // population description
                       field_strings[8],  // super_population
                       field_strings[9],  // super_description
                       field_strings[10],  // relationship
                       field_strings[11],  // siblings
                       field_strings[12],  // second_order
                       field_strings[13], // third_order
                       field_strings[14]); // comments

  if (not ped_data_->addPEDRecord(ped_record)) {

    ExecEnv::log().error("ParsePedFile::moveToPEDRecord; PED file: {}, record: {} cannot add PED record to map", fileName(), line_record);
    return false;

  }

  return true;

}

