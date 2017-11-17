//
// Created by kellerberrin on 17/11/17.
//



#include "kgl_alphabet_extend.h"


namespace kgl = kellerberrin::genome;


const kgl::ContigOffset_t kgl::ExtendDNA5::NUCLEOTIDE_COLUMNS;

const kgl::Nucleotide_t kgl::ExtendDNA5::A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendDNA5::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendDNA5::N_NUCLEOTIDE;
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
      ExecEnv::log().error("ExtendDNA5::Alphabet(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}

