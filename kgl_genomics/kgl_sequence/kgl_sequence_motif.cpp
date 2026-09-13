//
// Created by kellerberrin on 11/12/23.
//

#include "kgl_sequence_motif.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


namespace kellerberrin::genome::detail {   //  IUPAC code to regex string mapping.


// The code is mapped to its regex fragment; unmapped codes fall through to the default warn (identical output).
inline constexpr std::array<std::pair<char, std::string_view>, 17> IUPAC_REGEX_TABLE {{
  { 'A', "A" },    // Adenine
  { 'C', "C" },    // Cytosine
  { 'G', "G" },    // Guanine
  { 'T', "[TU]" }, // Thymine (or Uracil)
  { 'U', "[TU]" }, // Thymine (or Uracil)
  { 'R', "[AG]" }, // A or G
  { 'Y', "[CT]" }, // C or T
  { 'S', "[GC]" }, // G or C
  { 'W', "[AT]" }, // A or T
  { 'K', "[GT]" }, // G or T
  { 'M', "[AC]" }, // A or C
  { 'B', "[CGT]"}, // C or G or T
  { 'D', "[AGT]"}, // A or G or T
  { 'H', "[ACT]"}, // A or C or T
  { 'V', "[ACG]"}, // A or C or G
  { 'N', "." },    // any
  { '-', "." }     // any
}};


}   // end namespace


std::string kgl::SearchSequence::IUPACRegex(std::string_view IUPAC_search) {

  std::string upper_search = Utility::toupper(std::string(IUPAC_search));
  std::string regex_str;
  regex_str.reserve(upper_search.size() * 3);

  for (const char c : upper_search) {

    const auto table_iter = std::ranges::find(detail::IUPAC_REGEX_TABLE, c, &std::pair<char, std::string_view>::first);
    if (table_iter != detail::IUPAC_REGEX_TABLE.end()) {

      regex_str += table_iter->second;

    } else if (c == '.') {

      // Missing or any - the base is optional in the regex.
      regex_str += "?";

    } else {

      ExecEnv::log().warn("Non IUPAC nucleotide code: {} encountered in search string: {} - ignored", c, upper_search);

    }

  }

  return regex_str;

}