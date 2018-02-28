//
// Created by kellerberrin on 28/02/18.
//


#include "kgl_variant_factory_vcf_parse_impl.h"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


bool kgl::ParseVCFMiscImpl::parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                           const seqan::VcfHeader& header,
                                           ActiveContigMap& active_contig_map,
                                           bool cigar_required) {

  active_contig_map.clear();
  bool has_cigar = false;

  for (size_t idx = 0; idx != seqan::length(header); ++idx) {

    std::string key = seqan::toCString(header[idx].key);
    std::transform(key.begin(), key.end(), key.begin(), ::toupper);

    if (key == HEADER_CONTIG_KEY_) {

      std::map<std::string, std::string> item_map;
      std::string value = seqan::toCString(header[idx].value);
      tokenizeVcfHeaderKeyValues(value, item_map);

      auto id_result = item_map.find(ID_KEY_);
      auto length_result = item_map.find(CONTIG_LENGTH_KEY_);

      if (id_result != item_map.end() and length_result != item_map.end()) {

        ContigId_t contig_id = id_result->second;
        ContigSize_t contig_size = std::atoll(length_result->second.c_str());

        std::shared_ptr<const ContigFeatures> contig_ptr;
        if (genome_db_ptr->getContigSequence(contig_id, contig_ptr)) {

          if (contig_ptr->sequence().length() == contig_size) {

            active_contig_map[contig_id] = contig_size;

          } else {

            ExecEnv::log().warn("VCF header. VCF contig: {} size: {} not equal to the Genome database contig size: {}",
                                contig_id, contig_size, contig_ptr->sequence().length());
            ExecEnv::log().warn("Check that the VCF fasta file is compatible with the Genome database fasta file.");
          }

        } else {

          ExecEnv::log().warn("VCF header. VCF contig: {} size: {} not found in the genome database. ",
                              contig_id, contig_size);

        }

      }

    } else if (key == HEADER_INFO_KEY_) {

      std::map<std::string, std::string> item_map;
      std::string value = seqan::toCString(header[idx].value);
      tokenizeVcfHeaderKeyValues(value, item_map);

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
bool kgl::ParseVCFMiscImpl::tokenizeVcfHeaderKeyValues(const std::string& key_value_text,
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

      item_vec.push_back(seqan::toCString(*iter_item));

    }

    if (item_vec.size() >= 2) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = item_vec[1];

    } else {

      ExecEnv::log().warn("Problem parsing item: {} in VCF header line: {}, expected <key=value,..> pairs", *iter, key_value_text);
      return false;

    }

  }

  return true;

}


bool kgl::ParseVCFMiscImpl::parseCigar(const std::string& cigar,
                                       size_t& check_reference_size,
                                       size_t& check_alternate_size,
                                       std::vector<std::pair<char, size_t>>& parsed_cigar) {

  parsed_cigar.clear();
  check_reference_size = 0;
  check_alternate_size = 0;
  auto it = cigar.begin();
  char cigar_code;
  std::string cigar_size;
  while(it != cigar.end()) {

    cigar_size.clear();
    while(isdigit(*it) and it != cigar.end()) {

      cigar_size += *it;
      ++it;

    }

    if (cigar_size.empty()) {

      ExecEnv::log().error("VCF factory; cigar: {} contains no cigar size", cigar);
      return false;

    }

    size_t size = std::atoll(cigar_size.c_str());

    if (it != cigar.end()) {

      switch (*it) {

        case 'I':
          check_alternate_size += size;
          cigar_code = *it;
          break;

        case 'D':
          check_reference_size += size;
          cigar_code = *it;
          break;

        case 'X':
        case 'M':
          check_alternate_size += size;
          check_reference_size += size;
          cigar_code = *it;
          break;

        default:
          ExecEnv::log().error("VCF factory; cigar: {} contains unexpected cigar code: {}", cigar, *it);
          return false;

      }

      parsed_cigar.push_back(std::pair<char, size_t>(cigar_code, size));
      ++it;

    } else {

      ExecEnv::log().error("VCF factory; cigar: {} unexpected end; no terminating cigar code.", cigar);
      return false;

    }

  }

  return true;

}


// assumes input "key_1=value_1; ...;key_n=value_n"
bool kgl::ParseVCFMiscImpl::tokenizeVcfInfoKeyValues(const std::string& key_value_text,
                                                     std::map<std::string, std::string>& key_value_map) {

  key_value_map.clear();

  bt::char_separator<char> item_key_sep(";");
  bt::tokenizer<bt::char_separator<char>> tokenize(key_value_text, item_key_sep);
  for(auto iter = tokenize.begin(); iter != tokenize.end(); ++iter) {

    std::vector<std::string> item_vec;
    bt::char_separator<char> item_key_sep("=");
    bt::tokenizer<bt::char_separator<char>> tokenize_item(*iter, item_key_sep);
    for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

      item_vec.push_back(seqan::toCString(*iter_item));

    }

    if (item_vec.size() >= 2) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = item_vec[1];

    } else {

      ExecEnv::log().warn("Problem parsing item: {} in VCF Record line: {}, expected 'key=value;...;key=value' pairs", *iter, key_value_text);
      return false;

    }

  }

  return true;

}

