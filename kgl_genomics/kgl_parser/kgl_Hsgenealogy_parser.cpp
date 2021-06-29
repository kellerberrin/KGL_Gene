//
// Created by kellerberrin on 2/8/20.
//

#include "kgl_Hsgenealogy_parser.h"
#include "kel_utility.h"

#include <set>

namespace kgl = kellerberrin::genome;


bool kgl::HsGenomeGenealogyData::addGenealogyRecord(const HsGenealogyRecord& record) {

  auto [iter, result] = genealogy_record_map_.try_emplace(record.individualId(), record);

  if (not result) {

    ExecEnv::log().error("HsGenomeGenealogyData::addGenealogyRecord, could add PED for sample: {} (duplicate)", record.individualId());
    return false;

  }

  return true;

}


std::vector<kgl::GenomeId_t> kgl::HsGenomeGenealogyData::getGenomeList() const {

  std::vector<GenomeId_t> genome_list;
  for (auto const& [genome, genome_record] : getMap()) {

    genome_list.push_back(genome);

  }

  return genome_list;

}


std::optional<kgl::HsGenealogyRecord> kgl::HsGenomeGenealogyData::getGenomeGenealogyRecord(const std::string& genome) const {

  auto result = getMap().find(genome);
  if (result == getMap().end()) {

    return std::nullopt;

  }

  auto const& [genome_id, genome_record] = *result;

  return genome_record;

}


std::optional<kgl::HsGenomeAuxRecord> kgl::HsGenomeGenealogyData::getGenome(const std::string& genome) const {

  auto result = getMap().find(genome);
  if (result == getMap().end()) {

    return std::nullopt;

  }

  auto const& [genome_id, genome_record] = *result;

  return HsGenomeAuxRecord(genome_record);

}


// Create a list of populations.
void kgl::HsGenomeGenealogyData::refreshPopulationLists() {

  population_list_.clear();
  super_population_list_.clear();
  for (auto const& [genome_id, ped_record] : genealogy_record_map_) {

    population_list_[ped_record.population()] = ped_record.populationDescription();
    super_population_list_[ped_record.superPopulation()] = ped_record.superDescription();

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

void kgl::ParseHsGenomeGenealogyFile::readParseImpl(const std::string& file_name) {

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
  genealogy_data_->refreshPopulationLists();

}



bool kgl::ParseHsGenomeGenealogyFile::moveToRecord(const std::vector<std::string>& field_strings) {


  if (field_strings.size() != HsGenealogyRecord::genealogyFieldCount()) {

    ExecEnv::log().error("ParseHsGenomeGenealogyFile::moveToRecord; field count: {} not equal mandatory count: {}",
                         field_strings.size(), HsGenealogyRecord::genealogyFieldCount());
    return false;

  }

  HsGenealogyRecord ped_record(field_strings[0],   // family_id
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

  if (not genealogy_data_->addGenealogyRecord(ped_record)) {

    ExecEnv::log().error("ParseHsGenomeGenealogyFile::moveToRecord; record cannot add PED record to map");
    return false;

  }

  return true;

}

