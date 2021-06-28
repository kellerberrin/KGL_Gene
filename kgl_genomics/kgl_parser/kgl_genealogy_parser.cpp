//
// Created by kellerberrin on 2/8/20.
//

#include "kgl_genealogy_parser.h"
#include "kel_utility.h"

#include <set>

namespace kgl = kellerberrin::genome;


kgl::GenomeAuxRecord::GenomeAuxRecord(const GenealogyRecord& ped_record) {

  individual_id_ = ped_record.individualId();
  sex_ = ped_record.sexType();
  population_ = ped_record.population();
  population_description_ = ped_record.populationDescription();
  super_population_ = ped_record.superPopulation();
  super_description_ = ped_record.superDescription();

}




bool kgl::GenomeGenealogyData::addPEDRecord(const GenealogyRecord& record) {

  auto [iter, result] = genealogy_record_map_.try_emplace(record.individualId(), record);

  if (not result) {

    ExecEnv::log().error("GenomeGenealogyData::addPEDRecord, could add PED for sample: {} (duplicate)", record.individualId());
    return false;

  }

  return true;

}


std::vector<kgl::GenomeId_t> kgl::GenomeGenealogyData::getGenomeList() const {

  std::vector<GenomeId_t> genome_list;
  for (auto const& [genome, genome_record] : getMap()) {

    genome_list.push_back(genome);

  }

  return genome_list;

}


std::optional<kgl::GenealogyRecord> kgl::GenomeGenealogyData::getGenomePedRecord(const std::string& genome) const {

  auto result = getMap().find(genome);
  if (result == getMap().end()) {

    return std::nullopt;

  }

  auto const& [genome_id, genome_record] = *result;

  return genome_record;

}


std::optional<kgl::GenomeAuxRecord> kgl::GenomeGenealogyData::getGenome(const std::string& genome) const {

  auto result = getMap().find(genome);
  if (result == getMap().end()) {

    return std::nullopt;

  }

  auto const& [genome_id, genome_record] = *result;

  return GenomeAuxRecord(genome_record);

}


// Create a list of populations.
void kgl::GenomeGenealogyData::refreshPopulationLists() {

  population_list_.clear();
  super_population_list_.clear();
  for (auto const& [genome_id, ped_record] : genealogy_record_map_) {

    population_list_[ped_record.population()] = ped_record.populationDescription();
    super_population_list_[ped_record.superPopulation()] = ped_record.superDescription();

  }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

void kgl::ParseGenomeGenealogyFile::readParseImpl(const std::string& file_name) {

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



bool kgl::ParseGenomeGenealogyFile::moveToRecord(const std::vector<std::string>& field_strings) {


  if (field_strings.size() != GenealogyRecord::PEDFieldCount()) {

    ExecEnv::log().error("ParseGenomeGenealogyFile::moveToRecord; field count: {} not equal mandatory count: {}",
                         field_strings.size(), GenealogyRecord::PEDFieldCount());
    return false;

  }

  GenealogyRecord ped_record(field_strings[0],   // family_id
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

  if (not genealogy_data_->addPEDRecord(ped_record)) {

    ExecEnv::log().error("ParseGenomeGenealogyFile::moveToRecord; record cannot add PED record to map");
    return false;

  }

  return true;

}

