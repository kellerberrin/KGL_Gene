//
// Created by kellerberrin on 19/05/18.
//

#include "kgd_alleleread.h"
#include "kgd_utility.h"

#include <fstream>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

namespace bt = boost;
namespace kgd = kellerberrin::deconvolv;


bool kgd::AlleleReader::parseFile(const std::string& file_name, const std::string& delimiters) {

  std::ifstream data_file;
  size_t counter = 0;

  // Open input file.

  data_file.open(file_name);

  if (not data_file.good()) {

    ExecEnv::log().critical("AlleleReader; I/O error; could not open data file: {}", file_name);

  }

  try {

    std::string record_str;
    bt::char_separator<char> item_key_sep(delimiters.c_str());
    bool header_flag = false;

    while (not data_file.eof()) {

      std::vector<std::string> text_fields;

      std::getline(data_file, record_str);

      if(record_str.size() == 0) continue; // skip empty lines.

      if (record_str[0] == COMMENT_CHARACTER_) continue;   // ignore comment records.

      if (record_str[0] == HEADER_CHARACTER_) {

        header_flag = true;
        record_str.erase(0,1); // remove the header character.

      }

      bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);

      for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

        text_fields.push_back(Utility::trimWhiteSpace(*iter_item));

      }

      if (header_flag) {

        if (not parseHeaderLine(text_fields)) {

          ExecEnv::log().error("AlleleReader, Problem parsing header, file: {}, line record: {}", file_name, record_str);

        }

        header_flag = false; // reset the header flag.

      } else {

        if (not parseDataLine(text_fields)) {

          ExecEnv::log().error("AlleleReader, Problem parsing data line, file: {}, line record: {}", file_name, record_str);

        }

      }

      ++counter;

    }

    data_file.close();

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("AlleleReader, File: {}, unexpected I/O exception: {}", file_name, e.what());

  }

  ExecEnv::log().info("AlleleReader, Parsed: {} lines from file: {}", counter, file_name);

  return true;

}



bool kgd::AlleleReader::parseHeaderLine(const std::vector<std::string>& text_fields) {

  bool result = true;

  if (text_fields.size() < MINIMUM_FIELD_COUNT_) {

    ExecEnv::log().error("AlleleReader, Header items below minimum count: {} (#chrom, offset, value", MINIMUM_FIELD_COUNT_);
    result = false;

  } else {

  // The value header items are assumed to be named genomes (possibly combined as in "plaf"), create these in an array.

    for (size_t idx = MINIMUM_FIELD_COUNT_-1; idx < text_fields.size(); ++idx) {

      std::string genome_name = text_fields[idx];

      genome_allele_vector_.push_back(GenomeAlleles(genome_name));

    }

  }

  headers_ = text_fields;

  return result;

}


bool kgd::AlleleReader::parseDataLine(const std::vector<std::string>& text_fields) {

  bool result = true;

  if (text_fields.size() != headers_.size()) {

    ExecEnv::log().error("AlleleReader, Mismatch between header size: {} and data line size: {}", headers_.size(), text_fields.size());
    return false;

  }

  try {

    auto field_iter = text_fields.begin();

    ContigId_t contig_id = *field_iter;

    ++field_iter;

    ContigOffset_t offset = std::stoull(*field_iter);

    ++field_iter;

    size_t array_offset = 0;

    while (field_iter != text_fields.end()) {

      AlleleFreq_t allele_freq = std::stod(*field_iter);

      if (array_offset >= genome_allele_vector_.size()) {

        ExecEnv::log().error("AlleleReader, Mismatch between data line size: {} and genome array size: {}",
                             text_fields.size(), genome_allele_vector_.size());
        return false;

      }

      genome_allele_vector_[array_offset].addAlleleFreq(contig_id, offset, allele_freq);


      ++field_iter;
      ++array_offset;

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().error("AlleleReader, Conversion problem: {}", e.what());
    result = false;

  }

  return result;

}

