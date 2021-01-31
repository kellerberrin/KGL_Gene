//
// Created by kellerberrin on 2/8/20.
//

#include "kgl_ped_parser.h"
#include "kel_utility.h"

#include <set>

namespace kgl = kellerberrin::genome;



bool kgl::GenomePEDData::addPEDRecord(const PEDRecord& record) {

  auto [iter, result] = PED_record_map_.try_emplace(record.individualId(), record);

  if (not result) {

    ExecEnv::log().error("GenomePEDData::addPEDRecord, could add PED for sample: {} (duplicate)", record.individualId());
    return false;

  }

  return true;

}



// Create a list of populations.
void kgl::GenomePEDData::refreshPopulationLists() {

  population_list_.clear();
  super_population_list_.clear();
  for (auto const& [genome_id, ped_record] : PED_record_map_) {

    population_list_[ped_record.population()] = ped_record.populationDescription();
    super_population_list_[ped_record.superPopulation()] = ped_record.superDescription();

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

void kgl::ParsePedFile::readParsePEDImpl(const std::string& file_name) {

  std::shared_ptr<SquareTextRows> parsed_text = text_parser_.parseFlatFile(file_name, PED_FIELD_DELIMITER_CHAR_);

  bool skip_header = true;
  for (auto const& record_fields : parsed_text->getRowVector()) {

    if (skip_header) {

      skip_header = false;

    } else {

      moveToPEDRecord(record_fields);

    }

  }

  // Can be re-refreshed at any time.
  ped_data_->refreshPopulationLists();

}



bool kgl::ParsePedFile::moveToPEDRecord(const std::vector<std::string>& field_strings) {


  if (field_strings.size() != PEDRecord::PEDFieldCount()) {

    ExecEnv::log().error("ParsePedFile::moveToPEDRecord; field count: {} not equal mandatory count: {}",
                         field_strings.size(), PEDRecord::PEDFieldCount());
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
                       field_strings[8],  // super_population_
                       field_strings[9],  // super_description
                       field_strings[10],  // relationship
                       field_strings[11],  // siblings
                       field_strings[12],  // second_order
                       field_strings[13], // third_order
                       field_strings[14]); // comments

  if (not ped_data_->addPEDRecord(ped_record)) {

    ExecEnv::log().error("ParsePedFile::moveToPEDRecord; record cannot add PED record to map");
    return false;

  }

  return true;

}

