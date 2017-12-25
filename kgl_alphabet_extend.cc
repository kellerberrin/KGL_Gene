//
// Created by kellerberrin on 25/12/17.
//




#include "kgl_alphabet_extend.h"


namespace kgl = kellerberrin::genome;



const kgl::ContigOffset_t kgl::ExtendCountColumns::NUCLEOTIDE_COLUMNS;

const kgl::Nucleotide_t kgl::ExtendCountColumns::A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendCountColumns::C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendCountColumns::G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendCountColumns::U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendCountColumns::T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendCountColumns::N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::N_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::ExtendCountColumns::DELETE_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::ExtendCountColumns::DELETE_NUCLEOTIDE_OFFSET;


// Covert from char to alphabet
kgl::ExtendCountColumns::Alphabet kgl::ExtendCountColumns::convertChar(char chr_base) {

  // Translate the nucleotide to an array column
  switch (std::toupper(chr_base)) {

    case A_NUCLEOTIDE: return Alphabet::A;

    case C_NUCLEOTIDE: return Alphabet::C;

    case G_NUCLEOTIDE: return Alphabet::G;

    case U_NUCLEOTIDE:
    case T_NUCLEOTIDE: return Alphabet::T;

    case DELETE_NUCLEOTIDE: return Alphabet::X;

    default:
      ExecEnv::log().error("ExtendCountColumns::Alphabet(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}

