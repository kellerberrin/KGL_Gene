//
// Created by kellerberrin on 29/6/21.
//

#include "kgl_Hsgenealogy_parser.h"


namespace kgl = kellerberrin::genome;


kgl::HsGenomeAuxRecord::HsGenomeAuxRecord(const HsGenealogyRecord& ped_record) {

  individual_id_ = ped_record.individualId();
  sex_ = ped_record.sexType();
  population_ = ped_record.population();
  population_description_ = ped_record.populationDescription();
  super_population_ = ped_record.superPopulation();
  super_description_ = ped_record.superDescription();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information implementation object.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::HsGenomeAuxData::addGenomeAuxRecord(const HsGenomeAuxRecord& record) {

  auto [iter, result] = genome_aux_record_map_.try_emplace(record.individualId(), record);

  if (not result) {

    ExecEnv::log().error("HHsGenomeAuxData::addGenomeAuxRecord, could not add HsGenomeAuxRecord for sample: {} (duplicate)", record.individualId());
    return false;

  }

  return true;

}


std::vector<kgl::GenomeId_t> kgl::HsGenomeAuxData::getGenomeList() const {

  std::vector<GenomeId_t> genome_list;
  for (auto const& [genome, genome_record] : getMap()) {

    genome_list.push_back(genome);

  }

  return genome_list;

}


std::optional<kgl::HsGenomeAuxRecord> kgl::HsGenomeAuxData::getGenome(const std::string& genome) const {

  auto result = getMap().find(genome);
  if (result == getMap().end()) {

    return std::nullopt;

  }

  auto const& [genome_id, genome_record] = *result;

  return HsGenomeAuxRecord(genome_record);

}


// Create a list of populations.
void kgl::HsGenomeAuxData::refreshPopulationLists() {

  population_list_.clear();
  super_population_list_.clear();
  for (auto const& [genome_id, genome_aux_record] : genome_aux_record_map_) {

    population_list_[genome_aux_record.population()] = genome_aux_record.populationDescription();
    super_population_list_[genome_aux_record.superPopulation()] = genome_aux_record.superDescription();

  }

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// The Genome information parser object.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::ParseHsGenomeAuxFile::readParseImpl(const std::string& file_name) {

  std::shared_ptr<SquareTextRows> parsed_text = text_parser_.parseFlatFile(file_name, FIELD_DELIMITER_CHAR_);

  bool skip_header = true;
  for (auto const& record_fields : parsed_text->getRowVector()) {

    if (skip_header) {

      skip_header = false;

    } else {

      moveToRecord(record_fields);

    }

  }

  // Can be re-refreshed at any time.
  genome_aux_data_->refreshPopulationLists();

}



bool kgl::ParseHsGenomeAuxFile::moveToRecord(const std::vector<std::string>& field_strings) {


  if (field_strings.size() != GENOME_AUX_FIELD_COUNT_) {

    ExecEnv::log().error("ParseHsGenomeAuxFile::moveToRecord; field count: {} not equal mandatory count: {}",
                         field_strings.size(), GENOME_AUX_FIELD_COUNT_);
    return false;

  }

  HsGenomeAuxRecord genome_aux_record( field_strings[0],  //individual_id
                                       field_strings[1],  // sex
                                       field_strings[3],  // super_population_
                                       field_strings[4],  // super_description
                                       field_strings[5],  // population
                                       field_strings[6]);  // population description

  if (not genome_aux_data_->addGenomeAuxRecord(genome_aux_record)) {

    ExecEnv::log().error("ParseHsGenomeAuxFile::moveToRecord; record cannot add Genome Aux record to map");
    return false;

  }

  return true;

}

