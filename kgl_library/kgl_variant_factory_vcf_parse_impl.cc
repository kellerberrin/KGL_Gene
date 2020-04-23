//
// Created by kellerberrin on 28/02/18.
//

#include "kel_utility.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <edlib.h>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


namespace kgl = kellerberrin::genome;
namespace bt = boost;


bool kgl::ParseVCFMiscImpl::parseVcfHeader(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
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

      item_vec.push_back(*iter_item);

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
                                       std::vector<CigarEditItem>& parsed_cigar) {

  parsed_cigar.clear();
  check_reference_size = 0;
  check_alternate_size = 0;
  auto it = cigar.begin();
  CigarEditType cigar_code;
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
          cigar_code = CigarEditType::INSERT;
          break;

        case 'D':
          check_reference_size += size;
          cigar_code = CigarEditType::DELETE;
          break;

        case 'X':
          check_alternate_size += size;
          check_reference_size += size;
          cigar_code = CigarEditType::CHANGED;
          break;

        case 'M':
          check_alternate_size += size;
          check_reference_size += size;
          cigar_code = CigarEditType::UNCHANGED;
          break;

        default:
          ExecEnv::log().error("VCF factory; cigar: {} contains unexpected cigar code: {}", cigar, *it);
          return false;

      }

      parsed_cigar.push_back(CigarEditItem(size, cigar_code));
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

      item_vec.push_back(*iter_item);

    }

    if (item_vec.size() >= 2) {

      std::string item_key = item_vec[0];
      std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
      key_value_map[item_key] = item_vec[1];

    } else {

      if (item_vec.size() == 1) {

        std::string item_key = item_vec[0];
        std::transform(item_key.begin(), item_key.end(), item_key.begin(), ::toupper);
        key_value_map[item_key] = "";

      } else {

        ExecEnv::log().warn("Problem parsing item: {} in VCF Record line: {}, expected 'key=value;...;key=value' pairs", *iter, key_value_text);
        return false;

      }

    }

  }

  return true;

}


// assumes input "key_1=value_1; ...;key_n=value_n"
bool kgl::ParseVCFMiscImpl::tokenize(const std::string& parse_text,
                                     const std::string& separator_text,
                                     std::vector<std::string>& item_vector) {

  item_vector.clear();
  bt::char_separator<char> item_seperator(separator_text.c_str());
  bt::tokenizer<bt::char_separator<char>> tokenize_item(parse_text, item_seperator);
  for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

    item_vector.push_back(*iter_item);

  }

  return true;

}

// assumes input "key_1=value_1; ...;key_n=value_n"
bool kgl::ParseVCFMiscImpl::tokenize(std::string&& parse_text,
                                     const std::string& separator_text,
                                     std::vector<std::string>& item_vector) {

  item_vector.clear();
  bt::char_separator<char> item_seperator(separator_text.c_str());
  bt::tokenizer<bt::char_separator<char>> tokenize_item(parse_text, item_seperator);
  for(auto iter_item : tokenize_item) {

    item_vector.push_back(std::move(iter_item));

  }

  return true;

}

// Use edlib to generate a cigar string.
std::string kgl::ParseVCFMiscImpl::generateCigar(const CigarVector& cigar_vector) {

  std::stringstream ss;
  for(auto item : cigar_vector) {

    ss << item.first << static_cast<char>(item.second);

  }

  return ss.str();

}


// Use edlib to generate a cigar string.
std::string kgl::ParseVCFMiscImpl::generateCigar(const std::string& reference, const std::string& alternate) {

  return generateCigar(generateEditVector(reference, alternate));

}


// Use edlib to generate a cigar vector.
std::vector<kgl::CigarEditItem> kgl::ParseVCFMiscImpl::generateEditVector(const std::string& reference,
                                                                          const std::string& alternate) {

  std::vector<kgl::CigarEditItem> item_vector;
  std::vector<CigarEditType> edit_string;
  generateEditString(reference, alternate, edit_string);

  size_t same_count = 0;
  bool first_pass = true;
  CigarEditType previous_edit_item;
  for(auto edit : edit_string) {

    if (first_pass) {

      first_pass = false;
      previous_edit_item = edit;

    }

    if (previous_edit_item == edit) {

      ++same_count;

    } else {

      item_vector.emplace_back(CigarEditItem(same_count, previous_edit_item));
      same_count = 1;
      previous_edit_item = edit;

    }

  }

  if (not first_pass) {

    item_vector.emplace_back(CigarEditItem(same_count, previous_edit_item));

  }

  return item_vector;

}



// Use edlib to generate a cigar string.
void kgl::ParseVCFMiscImpl::generateEditString(const std::string& reference,
                                                 const std::string& alternate,
                                                 std::vector<CigarEditType>& edit_vector) {


  edit_vector.clear();

  EdlibAlignResult result = edlibAlign(alternate.c_str(), alternate.size(),reference.c_str(), reference.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  if (result.status == EDLIB_STATUS_OK) {

    for (int i = 0; i < result.alignmentLength; ++i) {

      switch(result.alignment[i]) {

        case 0:
          edit_vector.push_back(CigarEditType::UNCHANGED);
          break;

        case 1:
          edit_vector.push_back(CigarEditType::INSERT);
          break;

        case 2:
          edit_vector.push_back(CigarEditType::DELETE);
          break;

        case 3:
          edit_vector.push_back(CigarEditType::CHANGED);
          break;

        default:
          ExecEnv::log().error("ParseVCFMiscImpl::generateEditString; Unknown Cigar Code: {}, reference: {}, alternate: {}",
                               result.alignment[i], reference, alternate);
          break;

      }


    }


  } else {

    ExecEnv::log().error("ParseVCFMiscImpl::generateEditString; problem generating cigar reference:{}, alternate: {}", reference, alternate);

  }

  edlibFreeAlignResult(result);

}


// Given a reference count and a cigar vector compute a number that calculates the equivalent
// size of the alternate string.
// For UNCHANGED = 'M' and CHANGED = 'X' cigar items the reference_count and alternate count are incremented.
// For INSERT = 'I' the alternate is incremented and the reference_count is not.
// For DELETE = 'D' the reference count is incremented and the alternate is not.
size_t kgl::ParseVCFMiscImpl::alternateCount(size_t reference_count, const CigarVector& cigar_vector) {

  size_t reference_counter = 0;
  size_t alternate_counter = 0;
  auto cigar_ptr = cigar_vector.begin();

  while (reference_counter < reference_count and cigar_ptr != cigar_vector.end()) {

    switch(cigar_ptr->second) {

      case CigarEditType::UNCHANGED:
      case CigarEditType::CHANGED:
        if ((reference_counter + cigar_ptr->first) > reference_count) {

          alternate_counter += (reference_count - reference_counter);
          reference_counter = reference_count;

        } else {

          reference_counter += cigar_ptr->first;
          alternate_counter += cigar_ptr->first;

        }
        break;

      case CigarEditType::INSERT:
        alternate_counter += cigar_ptr->first;
        break;

      case CigarEditType::DELETE:
        if ((reference_counter + cigar_ptr->first) > reference_count) {

          reference_counter = reference_count;

        } else {

          reference_counter += cigar_ptr->first;

        }
        break;

    }

    ++cigar_ptr;

  }

  if (reference_counter != reference_count) {

    ExecEnv::log().error("ParseVCFMiscImpl::alternateCount(), unable to calculate alternate count for reference count: {}, cigar: {}",
                         reference_count, generateCigar(cigar_vector));
    alternate_counter = 0;

  }

  return alternate_counter;

}
