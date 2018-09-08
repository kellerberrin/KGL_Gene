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
  for (auto item : lc_aux_data_header) {

    aux_data_header_.push_back(Utility::trimWhiteSpace(Utility::toupper(item)));

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


size_t kgl::GenomeAuxData::fieldOffset(const std::string& field_name) const {

  std::string uc_field_name = Utility::trimWhiteSpace(Utility::toupper(field_name));

  size_t field_offset = 0;
  for (auto header_field : aux_data_header_) {

    if (uc_field_name == header_field) {

      return field_offset;

    }

    ++field_offset;

  }

  ExecEnv::log().warn("GenomeAuxData; could not find field: {}", uc_field_name);

  return 0;

}


bool kgl::GenomeAuxData::getGenomeAttributes(const std::string& genome_name, AuxAttributeVector& attribute_vector) const {

  auto result = aux_sample_information_.find(genome_name);

  if (result == aux_sample_information_.end()) {

    ExecEnv::log().error("GenomeAuxData; Genome: {} not found in auxiliary data", genome_name);
    return false;

  }

  attribute_vector = result->second;

  return true;

}

std::string kgl::GenomeAuxData::locationDate(const std::string& genome_name) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    ExecEnv::log().error("Could not find genome: {} in auxiliary data", genome_name);
    return "GENOME NOT FOUND";

  }

  std::string Country = attribute_vector[fieldOffset(COUNTRY)];
  std::string Site = attribute_vector[fieldOffset(SITE)];
  std::string Year = attribute_vector[fieldOffset(COLLECTION_YEAR)];

  return Country + " " + Site + " " + Year;

}


bool kgl::GenomeAuxData::locationDate(const std::string& genome_name,
                                      std::string& country,
                                      std::string& location,
                                      std::string& year) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    return false;

  }

  country = attribute_vector[fieldOffset(COUNTRY)];
  location = attribute_vector[fieldOffset(SITE)];
  year = attribute_vector[fieldOffset(COLLECTION_YEAR)];

  return true;

}



bool kgl::GenomeAuxData::isFieldSample(const std::string& genome_name) const {

  AuxAttributeVector attribute_vector;
  if (not getGenomeAttributes(genome_name, attribute_vector)) {

    ExecEnv::log().error("Could not find genome: {} in auxiliary data", genome_name);
    return false;

  }

  return attribute_vector[fieldOffset(FIELD_SAMPLE)] == TRUE_FLAG;

}
