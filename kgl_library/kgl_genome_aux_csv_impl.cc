//
// Created by kellerberrin on 19/03/18.
//


#include <fstream>
#include "kgl_genome_aux_csv.h"
#include "kgl_exec_env.h"
#include "kgl_variant_factory_vcf_parse_impl.h"


#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


bool kgl::GenomeAuxData::readParseAuxData(const std::string& aux_file_name) {


  std::ifstream aux_file;
  size_t counter = 0;

  // Open input file.

  aux_file.open(aux_file_name);

  if (not aux_file.good()) {

    ExecEnv::log().critical("GenomeAuxData; I/O error; could not open auxillary file: {}", aux_file_name);

  }

  try {


    bool first_line = true;
    while (true) {

      std::string record_str;

      if (std::getline(aux_file, record_str).eof()) break;

      if (not record_str.empty()) {

        if (record_str[0] == COMMENT) {

          continue;  // Skip comment lines.

        }

      }

      if (first_line) {

        if (not parseHeader(record_str)) {

          ExecEnv::log().error("GenomeAuxData, Problem parsing header (first) record: {}", record_str);

        }

        first_line = false;
        continue;

      }

      if (not parseDataline(record_str)) {

        ExecEnv::log().error("GenomeAuxData, Problem parsing data line record: {}", record_str);

      }

      ++counter;

    }

    aux_file.close();

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("GenomeAuxData, AUX file: {}, unexpected I/O exception: {}", aux_file_name, e.what());

  }

  ExecEnv::log().info("GenomeAuxData, Parsed: {} lines from aux file: {}", counter, aux_file_name);

  return true;

}


bool kgl::GenomeAuxData::parseHeader(const std::string& record_str) {

  AuxAttributeVector lc_aux_data_header;
  tokenize( record_str, lc_aux_data_header);

  // Covert Header to trimmed uppercase.
  for (const auto& item : lc_aux_data_header) {

    aux_data_header_.emplace_back(Utility::trimWhiteSpace(Utility::toupper(item)));

  }

  return true;

}


bool kgl::GenomeAuxData::parseDataline(const std::string& record_str) {

  AuxAttributeVector attributes;
  tokenize( record_str, attributes);

  if (attributes.size() != aux_data_header_.size() or attributes.size() == 0) {

    ExecEnv::log().error("GenomeAuxData; Header size: {} Aux size: {} mismatch; record line: {}", aux_data_header_.size(), attributes.size(), record_str);
    return false;
  }

  auto result = aux_sample_information_.insert(std::pair<std::string, AuxAttributeVector>(attributes[0], attributes));

  if (not result.second) {

    ExecEnv::log().error("GenomeAuxData; Duplicate record for item: {}; record line: {}", attributes[0], record_str);
    return false;

  }

  return true;

}


bool kgl::GenomeAuxData::tokenize(const std::string& parse_text,
                                  AuxAttributeVector& attribute_vector) {

  attribute_vector.clear();
  typedef boost::tokenizer<boost::escaped_list_separator<char>> tokenizer;
  tokenizer tokenize_item{parse_text};
  for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

    attribute_vector.push_back(Utility::trimWhiteSpace(*iter_item));

  }

  return true;

}

// Returns the field offset or boolean false if field not found.
bool kgl::GenomeAuxData::fieldOffset(const std::string& field_name, size_t& field_offset) const {

  std::string uc_field_name = Utility::trimWhiteSpace(Utility::toupper(field_name));

  for (size_t idx = 0; idx < aux_data_header_.size(); ++idx) {

    if (uc_field_name == aux_data_header_[idx]) {

      field_offset = idx;
      return true;

    }

  }

  ExecEnv::log().warn("GenomeAuxData; could not find field: {}", uc_field_name);

  field_offset = FIELD_NOT_FOUND;
  return false;

}

std::string kgl::GenomeAuxData::attributeValue(AuxAttributeVector& attribute_vector, const std::string& attribute_name) const {

  size_t field_offset;

  if (fieldOffset(attribute_name, field_offset)) {

    if (field_offset >= attribute_vector.size()) {

      ExecEnv::log().error("GenomeAuxData::attributeValue; field offset: {} for attribute: {} exceeds attribute vector size: {}",
                           field_offset, attribute_name, attribute_vector.size());
      std::string Empty;
      return Empty;

    }

    return attribute_vector[field_offset];

  }

  std::string Empty;
  return Empty;

}



bool kgl::GenomeAuxData::getGenomeAttributes(const GenomeId_t& genome_name, AuxAttributeVector& attribute_vector) const {

  auto result = aux_sample_information_.find(genome_name);

  if (result == aux_sample_information_.end()) {

    ExecEnv::log().error("GenomeAuxData; Genome: {} not found in auxiliary data", genome_name);
    return false;

  }

  attribute_vector = result->second;

  return true;

}

std::string kgl::GenomeAuxData::locationDate(const GenomeId_t& genome_name) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    ExecEnv::log().error("Could not find genome: {} in auxiliary data", genome_name);
    return "GENOME NOT FOUND";

  }

  return attributeValue(attribute_vector, COUNTRY) + " "
         + attributeValue(attribute_vector, SITE) + " "
         + attributeValue(attribute_vector, COLLECTION_YEAR);

}


bool kgl::GenomeAuxData::locationDate(const GenomeId_t& genome_name,
                                      std::string& country,
                                      std::string& location,
                                      std::string& year) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    return false;

  }

  country = attributeValue(attribute_vector, COUNTRY);
  location = attributeValue(attribute_vector, SITE);
  year = attributeValue(attribute_vector, COLLECTION_YEAR);

  return true;

}



kgl::GenomeId_t kgl::GenomeAuxData::source(const GenomeId_t& genome_name) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    ExecEnv::log().error("Could not find genome: {} in auxiliary data", genome_name);
    GenomeId_t empty;
    return empty;

  }

  return attributeValue(attribute_vector, SOURCE_GENOME);

}


bool kgl::GenomeAuxData::isFieldSample(const GenomeId_t& genome_name) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    ExecEnv::log().error("Could not find genome: {} in auxiliary data", genome_name);
    return false;

  }

  return attributeValue(attribute_vector, FIELD_SAMPLE) == TRUE_FLAG;

}


bool kgl::GenomeAuxData::isPreferredSample(const GenomeId_t& genome_name) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    ExecEnv::log().error("Could not find genome: {} in auxiliary data", genome_name);
    return false;

  }

  return attributeValue(attribute_vector, PREFERRED_SAMPLE) == TRUE_FLAG;

}


std::vector<std::string> kgl::GenomeAuxData::countryList() const {

  std::set<std::string> unique_country;

  for (auto genomes : getMap()) {

    std::string country = attributeValue(genomes.second, COUNTRY);
    unique_country.insert(country);

  }

  std::vector<std::string> country_list;
  for (auto country : unique_country) {

    country_list.push_back(country);

  }

  return country_list;

}


std::vector<kgl::GenomeId_t> kgl::GenomeAuxData::getCountry(const std::string& country) const {

  std::vector<GenomeId_t> genome_list;

  for (auto genomes : getMap()) {

    if (country == attributeValue(genomes.second, COUNTRY)) {

      genome_list.push_back(genomes.first);

    }

  }

  return genome_list;

}


std::vector<kgl::GenomeId_t> kgl::GenomeAuxData::getFieldSamples() const {

  std::vector<GenomeId_t> genome_list;

  for (auto genomes : getMap()) {

    if (isFieldSample(genomes.first)) {

      genome_list.push_back(genomes.first);

    }

  }

  return genome_list;

}


std::vector<kgl::GenomeId_t> kgl::GenomeAuxData::getPreferredSamples() const {

  std::vector<GenomeId_t> genome_list;

  for (auto genomes : getMap()) {

    if (isPreferredSample(genomes.first)) {

      genome_list.push_back(genomes.first);

    }

  }

  return genome_list;

}


// Convenience static function splits the phased populations into different countries (preferred genomes only).
std::vector<kgl::CountryPair> kgl::GenomeAuxData::getCountries(const std::string& aux_file,
                                                               std::shared_ptr<const kgl::PhasedPopulation> population_ptr) {

  // Read the AUX data.
  kgl::GenomeAuxData aux_data;
  aux_data.readParseAuxData(aux_file);

  // Get a list of countries.
  std::vector<std::string> countries = aux_data.countryList();

  // Only preferred genomes.
  std::shared_ptr<const kgl::PhasedPopulation> preferred_pop_ptr = population_ptr->filterGenomes("Preferred", aux_data.getPreferredSamples());

  std::vector<CountryPair> pair_vector;
  // Create a vector of all countries
  for (auto country : countries) {

    std::vector<GenomeId_t> genome_list = aux_data.getCountry(country);
    // First element is the genome name; second element is the genome name in the (VCF) population.
    std::vector<std::pair<GenomeId_t, GenomeId_t >> source_pairs;

    for (auto genome : genome_list) {

      GenomeId_t vcf_source = aux_data.source(genome);

      if (not vcf_source.empty()) {

        source_pairs.emplace_back(std::pair<GenomeId_t, GenomeId_t >(genome, vcf_source));

      } else {

        ExecEnv::log().warn("getCountries; Could not find source (vcf genome id) for genome id: {}", genome);

      }

    }

    std::shared_ptr<const kgl::PhasedPopulation> country_pop_ptr = preferred_pop_ptr->filterRenameGenomes(country, source_pairs);

    // Ignore empty populations
    if (not country_pop_ptr->getMap().empty()) {

      pair_vector.emplace_back(CountryPair(country, country_pop_ptr));

    }

  }

  return pair_vector;

}
