//
// kgl_sequence_motif.cpp — IUPAC code conversion (consteval table).
//

#include "kgl_sequence_motif.h"
#include "kel_utility.h"
#include "kel_exec_env.h"

#include <array>
#include <cctype>


namespace kgl = kellerberrin::genome;


namespace kellerberrin::genome::detail {   //  IUPAC code to regex string mapping.


// The code is mapped to its regex fragment; unmapped codes fall through to the default warn
// (identical output to the reference). Indexed by the upper-cased raw char.
inline constexpr std::array<std::string_view, 256> IUPAC_REGEX_TABLE = []() consteval {

  std::array<std::string_view, 256> table{};
  table.fill(std::string_view{});
  table[static_cast<unsigned char>('A')] = "A";
  table[static_cast<unsigned char>('C')] = "C";
  table[static_cast<unsigned char>('G')] = "G";
  table[static_cast<unsigned char>('T')] = "[TU]";
  table[static_cast<unsigned char>('U')] = "[TU]";
  table[static_cast<unsigned char>('R')] = "[AG]";
  table[static_cast<unsigned char>('Y')] = "[CT]";
  table[static_cast<unsigned char>('S')] = "[GC]";
  table[static_cast<unsigned char>('W')] = "[AT]";
  table[static_cast<unsigned char>('K')] = "[GT]";
  table[static_cast<unsigned char>('M')] = "[AC]";
  table[static_cast<unsigned char>('B')] = "[CGT]";
  table[static_cast<unsigned char>('D')] = "[AGT]";
  table[static_cast<unsigned char>('H')] = "[ACT]";
  table[static_cast<unsigned char>('V')] = "[ACG]";
  table[static_cast<unsigned char>('N')] = ".";
  table[static_cast<unsigned char>('-')] = ".";
  return table;

}();


}   // end namespace


std::string kgl::SearchSequence::IUPACRegex(std::string_view IUPAC_search) {

  const std::string upper_search = Utility::toupper(std::string(IUPAC_search));
  std::string regex_str;
  regex_str.reserve(upper_search.size() * 3);

  for (const char c : upper_search) {

    const std::string_view fragment = detail::IUPAC_REGEX_TABLE[static_cast<unsigned char>(c)];
    if (not fragment.empty()) {

      regex_str += fragment;

    } else if (c == '.') {

      // Missing or any - the base is optional in the regex.
      regex_str += "?";

    } else {

      ExecEnv::log().warn("Non IUPAC nucleotide code: {} encountered in search string: {} - ignored", c, upper_search);

    }

  }

  return regex_str;

}
