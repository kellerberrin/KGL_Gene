//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_alphabet_base.h"

namespace kgl = kellerberrin::genome;

// Instantiate the ExtendDNA5 constants for the edification of the linker.
const kgl::Nucleotide_t kgl::DNA5::A_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::DNA5::A_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::DNA5::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::C_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::DNA5::C_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::DNA5::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::G_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::DNA5::G_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::DNA5::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::U_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::DNA5::U_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::DNA5::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::T_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::DNA5::T_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::DNA5::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::N_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::DNA5::N_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::DNA5::N_NUCLEOTIDE_OFFSET;


// Convert char to Alphabet enum type.
kgl::DNA5::Alphabet kgl::DNA5::convertChar(char chr_base) {

  switch (std::toupper(chr_base)) {

    case 'A':return Alphabet::A;

    case 'C': return Alphabet::C;

    case 'G': return Alphabet::G;

    case 'U':
    case 'T': return Alphabet::T;

    case 'N': return Alphabet::N;

    default:
      ExecEnv::log().error("BaseDNA5::convertchar(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}



const kgl::ContigOffset_t kgl::ExtendDNA5::NUCLEOTIDE_COLUMNS;

const kgl::Nucleotide_t kgl::ExtendDNA5::A_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::ExtendDNA5::A_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::ExtendDNA5::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::C_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::ExtendDNA5::C_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::ExtendDNA5::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::G_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::ExtendDNA5::G_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::ExtendDNA5::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::U_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::ExtendDNA5::U_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::ExtendDNA5::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::T_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::ExtendDNA5::T_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::ExtendDNA5::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::N_NUCLEOTIDE;
const kgl::Nucleotide_t kgl::ExtendDNA5::N_NUCLEOTIDE_LC;
const kgl::ContigOffset_t kgl::ExtendDNA5::N_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::DELETE_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::DELETE_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::INSERT_A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::INSERT_A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::INSERT_C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::INSERT_C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::INSERT_G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::INSERT_G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::INSERT_U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::INSERT_U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::INSERT_T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::INSERT_T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::INSERT_N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::INSERT_N_NUCLEOTIDE_OFFSET;


// Covert from char to alphabet
kgl::ExtendDNA5::Alphabet kgl::ExtendDNA5::convertChar(char chr_base) {

  // Translate the nucleotide to an array column
  switch (chr_base) {

    case A_NUCLEOTIDE:
    case A_NUCLEOTIDE_LC: return Alphabet::A;

    case C_NUCLEOTIDE:
    case C_NUCLEOTIDE_LC: return Alphabet::C;

    case G_NUCLEOTIDE:
    case G_NUCLEOTIDE_LC: return Alphabet::G;

    case T_NUCLEOTIDE:
    case T_NUCLEOTIDE_LC: return Alphabet::T;

    case U_NUCLEOTIDE:
    case U_NUCLEOTIDE_LC: return Alphabet::T;

    case DELETE_NUCLEOTIDE: return Alphabet::X;

    case INSERT_A_NUCLEOTIDE: return Alphabet::E;

    case INSERT_C_NUCLEOTIDE: return Alphabet::F;

    case INSERT_G_NUCLEOTIDE: return Alphabet::I;

    case INSERT_T_NUCLEOTIDE: return Alphabet::J;

    case INSERT_N_NUCLEOTIDE: return Alphabet::K;

    default:
      ExecEnv::log().error("ExtendDNA5::Alphabet(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}
