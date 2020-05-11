//
// Created by kellerberrin on 11/5/20.
//

#include "kgl_variant_factory_vcf_parse_header.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::VCFParseHeader::parseHeader(const std::string& vcf_file_name, std::unique_ptr<BaseStreamIO>&& vcf_stream) {

  // Open input file.
  if (not vcf_stream->open(vcf_file_name)) {

    ExecEnv::log().critical("I/O error; could not open VCF file: {}", vcf_file_name);

  }

  try {

    long counter = 0;
    bool found_header = false;

    while (true) {

      IOLineRecord line_record = vcf_stream->readLine();
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


bool kgl::VCFParseHeader::parseVcfHeader( std::shared_ptr<const GenomeReference> genome_db_ptr,
                                          const VcfHeaderInfo& header_info,
                                          ActiveContigMap& active_contig_map,
                                          bool cigar_required) {

  active_contig_map.clear();
  bool has_cigar = false;

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
        ContigSize_t contig_size = std::atoll(length_result->second.c_str());

        std::optional<std::shared_ptr<const ContigReference>> contig_opt = genome_db_ptr->getContigSequence(contig_id);
        if (contig_opt) {

          if (contig_opt.value()->sequence().length() == contig_size) {

            active_contig_map[contig_id] = contig_size;

          } else {

            ExecEnv::log().warn("VCF header. VCF contig: {} size: {} not equal to the Genome database contig size: {}",
                                contig_id, contig_size, contig_opt.value()->sequence().length());
            ExecEnv::log().warn("Check that the VCF fasta file is compatible with the Genome database fasta file.");
          }

        } else {

          ExecEnv::log().info("VCF header. VCF contig: {} size: {} not found in the genome database. ",
                              contig_id, contig_size);

        }

      }

    } else if (key == HEADER_INFO_KEY_) {

      std::map<std::string, std::string> item_map;

      if (not tokenizeVcfHeaderKeyValues(value, item_map)) {

        ExecEnv::log().warn("VCF header. Cannot tokenize VCF file header string: {}", value);

      }

      auto id_result = item_map.find(ID_KEY_);

      if (id_result != item_map.end()) {

        std::string value = id_result->second;
        std::transform(value.begin(), value.end(), value.begin(), ::toupper);
        if (value == ID_CIGAR_VALUE_) {

          has_cigar = true;

        }

      }

    }

  }

  if (active_contig_map.size() == genome_db_ptr->contigCount()) {

    ExecEnv::log().info("Genome database contigs: {} and VCF contigs: {}; all contig sizes match.",
                        genome_db_ptr->contigCount(), active_contig_map.size());

  } else {

    ExecEnv::log().warn("VCF file contig count: {}. Genome database contig count: {}", active_contig_map.size(), genome_db_ptr->contigCount());
    ExecEnv::log().warn("No variants will be generated for the missing contigs. The missing contigs are:");

    for (auto contig : genome_db_ptr->getMap()) {

      auto result = active_contig_map.find(contig.first);

      if (result == active_contig_map.end()) {

        ExecEnv::log().warn("Contig: {} present in genome database, missing from VCF file.", contig.first);

      }

    }

  }

  if ((not has_cigar) and cigar_required) {

    ExecEnv::log().error("This VCF file does not define a 'CIGAR' field; this field is required to parse 'freebayes' VCF files.");
    ExecEnv::log().error("See the VCF format of the variant caller 'freebayes' for more information.");
    return false;

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

