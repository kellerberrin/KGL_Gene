//
// Created by kellerberrin on 11/5/20.

#include "kgl_variant_factory_vcf_parse_header.h"
#include "kel_basic_io.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::VCFParseHeader::parseHeader(const std::string& vcf_file_name) {

  // Open input file. Plain text or compressed.
  std::optional<std::unique_ptr<BaseStreamIO>> vcf_stream_opt = BaseStreamIO::getReaderStream(vcf_file_name);
  if (not vcf_stream_opt) {

    ExecEnv::log().critical("I/O error; could not open VCF file: {}", vcf_file_name);

  }

  try {

    long counter = 0;
    bool found_header = false;

    while (true) {

      IOLineRecord line_record = vcf_stream_opt.value()->readLine();
      if (not line_record) break;

      std::string& record_str = *line_record.value().second;

      size_t pos = record_str.find_first_of(KEY_SEPARATOR_);

      if (pos != std::string::npos) {

        std::string key = record_str.substr(0, pos);
        std::string value = record_str.substr(pos, std::string::npos);

        pos = value.find_first_of(KEY_SEPARATOR_);

        if (pos != std::string::npos) {

          value = value.erase(pos, pos + std::string(KEY_SEPARATOR_).length());

        }

        pos = key.find_first_of(KEY_PREFIX_);

        if (pos != std::string::npos) {

          key = key.erase(pos, pos + std::string(KEY_PREFIX_).length());

        }

        vcf_header_info_.emplace_back(key, value);

      }

      std::string line_prefix = record_str.substr(0, FIELD_NAME_FRAGMENT_LENGTH_);
      if (line_prefix == FIELD_NAME_FRAGMENT_) {

        found_header = true;
        std::vector<std::string> field_vector = Utility::tokenizer(record_str, "\t");
        size_t field_count = 0;
        for(auto const& field :field_vector) {

          if (field_count >= SKIP_FIELD_NAMES_) {

            vcf_genomes_.push_back(field);

          }

          ++field_count;

        }

        break; // #CHROM is the last field in the VCF header so stop processing.

      }

      ++counter;

    }

    if (not found_header) {

      ExecEnv::log().error("VCF Genome Names Not Found");

    } else {

      ExecEnv::log().info("{} Genomes in VCF Header {}, Header lines processed: {}", vcf_genomes_.size(), vcf_file_name, counter);

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCFParseHeader::parseHeader; VCF file: {}, unexpected I/O exception: {}", vcf_file_name, e.what());

  }

  return true;

}


bool kgl::VCFParseHeader::parseVcfHeader(const VcfHeaderInfo& header_info, VCFContigMap& vcf_contig_map, VCFInfoRecordMap& vcf_info_map) {

  vcf_contig_map.clear();

  for (auto const& [lower_key, value] : header_info) {

    std::string key = Utility::toupper(lower_key);

    if (key == HEADER_CONTIG_KEY_) {

      std::map<std::string, std::string> item_map;

      if (not tokenizeVcfHeaderKeyValues(value, item_map)) {

        ExecEnv::log().warn("VCF header. Cannot tokenize VCF file header string: {}", value);

      }

      auto id_result = item_map.find(ID_KEY_);
      auto length_result = item_map.find(CONTIG_LENGTH_KEY_);

      if (id_result != item_map.end() and length_result != item_map.end()) {

        ContigId_t contig_id = id_result->second;
        ContigSize_t contig_size = std::stol(length_result->second);

        vcf_contig_map[contig_id] = contig_size;

      }

    } else if (key == HEADER_INFO_KEY_) {

      std::map<std::string, std::string> item_map;

      if (not tokenizeVcfHeaderKeyValues(value, item_map)) {

        ExecEnv::log().warn("VCF header. Cannot tokenize VCF file header string: {}", value);

      }

      VCFInfoRecord info_record;
      info_record.ID = item_map[ID_KEY_];
      info_record.description = item_map[DESCRIPTION_KEY_];
      info_record.number = item_map[NUMBER_KEY_];
      info_record.type = item_map[TYPE_KEY_];

      vcf_info_map[info_record.ID] = info_record;

    }

  }

  return true;

}

// assumes input "<key_1=value_1, ...,key_n=value_n>"
bool kgl::VCFParseHeader::tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
                                                       std::map<std::string, std::string>& key_value_map) {

  key_value_map.clear();
  std::string copy_key_value_text = key_value_text;
  // Remove the angle backets from the value.
  bt::erase_all(copy_key_value_text, "<");
  bt::erase_all(copy_key_value_text, ">");


  bt::tokenizer<bt::escaped_list_separator<char>> tokenize(copy_key_value_text);
  for(auto iter = tokenize.begin(); iter != tokenize.end(); ++iter) {

    std::vector<std::string> item_vec;
    bt::char_separator<char> item_key_sep("=");
    bt::tokenizer<bt::char_separator<char>> tokenize_item(*iter, item_key_sep);
    for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

      item_vec.push_back(*iter_item);

    }

    if (item_vec.size() >= 2) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = item_vec[1];

    } else if (item_vec.size() == 1) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = "";

    } else {

      ExecEnv::log().warn("Problem parsing item: {} in VCF header line: {}, expected <key=value,..> pairs", *iter, key_value_text);
      return false;

    }

  }

  return true;

}

bool kgl::VCFParseHeader::checkVCFReferenceContigs( const VCFContigMap& vcf_contig_map,
                                                    std::shared_ptr<const GenomeReference> reference_genome) {

  bool contigs_found = true;
  for (auto const& [contig_id, contig_ptr] : reference_genome->getMap()) {

    auto result = vcf_contig_map.find(contig_id);
    if (result ==  vcf_contig_map.end()) {

      contigs_found = false;
      break;

    }

    if (result->second != contig_ptr->sequence().length()) {

      ExecEnv::log().warn("VCFParseHeader::checkVCFReferenceContigs, Genome: {}, Contig: {}, mismatch in VCF size: {} and Reference Contig size: {}",
                          reference_genome->genomeId(), contig_id, result->second, contig_ptr->sequence().length());
      contigs_found = false;

    }

  }

  return contigs_found;

}


bool kgl::VCFParseHeader::VCFContigAliasRemapping(const VCFContigMap& vcf_contig_map,
                                                  std::shared_ptr<const GenomeReference> reference_genome,
                                                  VCFContigAliasMap contig_alias_map) {

  // Reverse the VCF contig map.
  std::map<size_t, std::string> reverse_map;
  for (auto const& [contig_id, contig_size] : vcf_contig_map) {

    reverse_map[contig_size] = contig_id;

  }

  // Now lookup all the reference genome contigs.
  for (auto const& [contig_id, contig_ptr] : reference_genome->getMap()) {

    auto result = reverse_map.find(contig_ptr->sequence().length());
    if (result == reverse_map.end()) {

      contig_alias_map[contig_id] = contig_id;

    } else {

      contig_alias_map[result->second] = contig_id;

    }

  }

  if (contig_alias_map.size() != vcf_contig_map.size()) {

    ExecEnv::log().warn("VCFParseHeader::VCFContigAliasRemapping, Caution - VCF/References Alias size: {} not equal to VCF contig size: {}",
                        contig_alias_map.size(), vcf_contig_map.size());

  }

  return true;

}


