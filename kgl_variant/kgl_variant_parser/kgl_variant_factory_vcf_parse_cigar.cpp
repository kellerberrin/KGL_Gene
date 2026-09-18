//
// Created by kellerberrin on 28/02/18.
//

#include "kel_utility.h"
#include "kel_exec_env.h"
#include "kgl_variant_factory_vcf_parse_cigar.h"

#include <edlib.h>

#include <boost/algorithm/string.hpp>
#include <sstream>


namespace kgl = kellerberrin::genome;
namespace bt = boost;



bool kgl::ParseVCFCigar::parseCigar(const std::string& cigar,
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

// Use edlib to generate a cigar string.
std::string kgl::ParseVCFCigar::generateCigar(const CigarVector& cigar_vector) {

  std::stringstream ss;
  for(auto item : cigar_vector) {

    ss << item.first << static_cast<char>(item.second);

  }

  return ss.str();

}


// Use edlib to generate a cigar string.
std::string kgl::ParseVCFCigar::generateCigar(const std::string& reference, const std::string& alternate) {

  return generateCigar(generateEditVector(reference, alternate));

}

// Calculate the allele offset, the non-zero value of the first CigarVector item if it is 'UNCHANGED'
kgl::ContigOffset_t kgl::ParseVCFCigar::alleleEditOffset(const std::string& reference, const std::string& alternate) {

  std::vector<kgl::CigarEditItem> cigar =  generateEditVector(reference, alternate);

  if (cigar.empty()) {

    return 0;

  } else if (cigar.front().second == CigarEditType::UNCHANGED) {

    return cigar.front().first;

  } else {

    return 0;

  }

}


// Calculate the allele offset, the non-zero value of the first CigarVector item if it is 'UNCHANGED'
kgl::ContigOffset_t kgl::ParseVCFCigar::alleleOffset(const std::string& reference, const std::string& alternate) {

  const size_t common_size = std::min(reference.size(), alternate.size());

  size_t common_prefix{0};
  for (size_t index = 0; index < common_size; ++index) {

    if (reference[index] != alternate[index]) {

      return common_prefix;

    }

    ++common_prefix;

  }

  return common_prefix;

}

// The CigarVector contains (n x 'M') + (1 x 'X') and is, in effect, an SNP.
bool kgl::ParseVCFCigar::isSNP(const std::string& reference, const std::string& alternate) {

  // referance and alternate are different sizes then cannot be an SNP.
  if (reference. size() != alternate.size()) {

    return false;

  }

  std::vector<kgl::CigarEditType> edit_vector = generateEditString(reference, alternate);

  size_t x_count{0};
  for (auto item : edit_vector) {

    if (item == CigarEditType::CHANGED) {

      // More than 1 'X'
      if (x_count > 0) {

        return false;

      }

      ++x_count;

    } else if (item != CigarEditType::UNCHANGED) {

      return false;

    }

  }

  // 1'X' and n'M'
  return true;

}


// Use edlib to generate a cigar vector.
std::vector<kgl::CigarEditItem> kgl::ParseVCFCigar::generateEditVector(const std::string& reference,
                                                                       const std::string& alternate) {

  std::vector<kgl::CigarEditItem> item_vector;

  std::vector<CigarEditType> edit_string = generateEditString(reference, alternate);

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

      item_vector.emplace_back(same_count, previous_edit_item);
      same_count = 1;
      previous_edit_item = edit;

    }

  }

  if (not first_pass) {

    item_vector.emplace_back(same_count, previous_edit_item);

  }

  return item_vector;

}



// Use edlib to generate a cigar string.
std::vector<kgl::CigarEditType> kgl::ParseVCFCigar::generateEditString(const std::string& reference,
                                                                       const std::string& alternate) {


  std::vector<CigarEditType> edit_vector;

  EdlibAlignResult result = edlibAlign(alternate.c_str(),
                                       alternate.size(),
                                       reference.c_str(),
                                       reference.size(),
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
          ExecEnv::log().error("ParseVCFCigar::generateEditString; Unknown Cigar Code: {}, reference: {}, alternate: {}",
                               result.alignment[i], reference, alternate);
          break;

      }


    }


  } else {

    ExecEnv::log().error("ParseVCFCigar::generateEditString; problem generating cigar reference:{}, alternate: {}", reference, alternate);

  }

  edlibFreeAlignResult(result);

  return edit_vector;

}


// Given a reference count and a cigar vector compute a number that calculates the equivalent
// size of the alternate string.
// For UNCHANGED = 'M' and CHANGED = 'X' cigar items the reference_count and alternate count are incremented.
// For INSERT = 'I' the alternate is incremented and the reference_count is not.
// For DELETE = 'D' the reference count is incremented and the alternate is not.
size_t kgl::ParseVCFCigar::alternateCount(size_t reference_count, const CigarVector& cigar_vector) {

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

    ExecEnv::log().error("ParseVCFCigar::alternateCount(), unable to calculate alternate count for reference count: {}, cigar: {}",
                         reference_count, generateCigar(cigar_vector));
    alternate_counter = 0;

  }

  return alternate_counter;

}

