//
// Created by kellerberrin on 11/12/23.
//

#include "kgl_sequence_motif.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


std::string kgl::SearchSequence::IUPACRegex(const std::string_view& IUPAC_search) {

  std::string upper_search = Utility::toupper(std::string(IUPAC_search));
  std::string regex_str;

  for (const char c : upper_search) {

    switch(c) {

      case 'A': // Adenine
        regex_str += "A";
        break;

      case 'C':	// Cytosine
        regex_str += "C";
        break;

      case 'G': //	Guanine
        regex_str += "G";
        break;

      case 'T': //	Thymine (or Uracil)
      case 'U': //	Thymine (or Uracil)
        regex_str += "[TU]";
        break;

      case 'R':	// A or G
        regex_str += "[AG]";
        break;

      case 'Y': //	C or T
        regex_str += "[CT]";
        break;

      case 'S': //	G or C
        regex_str += "[GC]";
        break;

      case 'W':	// A or T
        regex_str += "[AT]";
        break;

      case 'K': // G or T
        regex_str += "[GT]";
        break;

      case 'M':	// A or C
        regex_str += "[AC]";
        break;

      case 'B': //C or G or T
        regex_str += "[CGT]";
        break;

      case 'D':	// A or G or T
        regex_str += "[AGT]";
        break;

      case 'H': //	A or C or T
        regex_str += "[ACT]";
        break;

      case 'V': //	A or C or G
        regex_str += "[ACG]";
        break;

      case 'N': // any
      case '-': // any
        regex_str += ".";
        break;

      case '.': // missing or any
        regex_str += "?";
        break;

      default:
        ExecEnv::log().warn("Non IUPAC nucleotide code: {} encountered in search string: {} - ignored", c, upper_search);
        break;

    }


  }

  return regex_str;

}

