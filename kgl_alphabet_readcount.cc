//
// Created by kellerberrin on 24/12/17.
//




#include "kgl_alphabet_readcount.h"


namespace kgl = kellerberrin::genome;


const kgl::ContigOffset_t kgl::ReadCountColumns::NUCLEOTIDE_COLUMNS;

const kgl::Nucleotide_t kgl::ReadCountColumns::A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::N_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::DELETE_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::DELETE_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::INSERT_A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::INSERT_A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::INSERT_C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::INSERT_C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::INSERT_G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::INSERT_G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::INSERT_U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::INSERT_U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::INSERT_T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::INSERT_T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ReadCountColumns::INSERT_N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ReadCountColumns::INSERT_N_NUCLEOTIDE_OFFSET;


// Covert from char to alphabet
kgl::ReadCountColumns::Alphabet kgl::ReadCountColumns::convertChar(char chr_base) {

  // Translate the nucleotide to an array column
  switch (std::toupper(chr_base)) {

    case A_NUCLEOTIDE: return Alphabet::A;

    case C_NUCLEOTIDE: return Alphabet::C;

    case G_NUCLEOTIDE: return Alphabet::G;

    case U_NUCLEOTIDE:
    case T_NUCLEOTIDE: return Alphabet::T;

    case DELETE_NUCLEOTIDE: return Alphabet::X;

    case INSERT_A_NUCLEOTIDE: return Alphabet::E;

    case INSERT_C_NUCLEOTIDE: return Alphabet::F;

    case INSERT_G_NUCLEOTIDE: return Alphabet::I;

    case INSERT_T_NUCLEOTIDE: return Alphabet::J;

    case INSERT_N_NUCLEOTIDE: return Alphabet::K;

    default:
      ExecEnv::log().error("ReadCountColumns::Alphabet(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}


