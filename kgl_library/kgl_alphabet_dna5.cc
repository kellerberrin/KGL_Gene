//
// Created by kellerberrin on 17/11/17.
//


#include "kgl_alphabet_dna5.h"


namespace kgl = kellerberrin::genome;


const kgl::ContigOffset_t kgl::DNA5::NUCLEOTIDE_COLUMNS;

const kgl::Nucleotide_t kgl::DNA5::A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::DNA5::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::DNA5::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::DNA5::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::DNA5::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::DNA5::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::DNA5::N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::DNA5::N_NUCLEOTIDE_OFFSET;


bool kgl::DNA5::validAlphabet(Alphabet nucleotide) {

  auto int_value = static_cast<size_t>(nucleotide);

// We DO NOT use a switch here.
// Because the switch assumes we can only have 5 base types (and a memory corrupted sequence may not).
  bool compare = int_value == static_cast<size_t>(Alphabet::A)
                 or int_value == static_cast<size_t>(Alphabet::C)
                 or int_value == static_cast<size_t>(Alphabet::G)
                 or int_value == static_cast<size_t>(Alphabet::T)
                 or int_value == static_cast<size_t>(Alphabet::N);

  return compare;

}



// Convert char to Alphabet enum type.
kgl::DNA5::Alphabet kgl::DNA5::convertChar(char chr_base) {

  switch (std::toupper(chr_base)) {

    case A_NUCLEOTIDE:return Alphabet::A;

    case C_NUCLEOTIDE: return Alphabet::C;

    case G_NUCLEOTIDE: return Alphabet::G;

    case U_NUCLEOTIDE:
    case T_NUCLEOTIDE: return Alphabet::T;

    case N_NUCLEOTIDE: return Alphabet::N;

    default:
      ExecEnv::log().error("DNA5::convertchar(), Invalid nucleotide: {}", chr_base);
      return Alphabet::N;

  }

}


// Find complementary bases.
kgl::CodingDNA5::Alphabet kgl::DNA5::complementNucleotide(Alphabet nucleotide) {

  // Translate the nucleotide
  switch (nucleotide) {

    case Alphabet::A: return CodingDNA5::Alphabet::T;
    case Alphabet::C: return CodingDNA5::Alphabet::G;
    case Alphabet::G: return CodingDNA5::Alphabet::C;
    case Alphabet::T: return CodingDNA5::Alphabet::A;
    case Alphabet::N: return CodingDNA5::Alphabet::N;

  }

  return CodingDNA5::Alphabet::N; //  Never reached, to keep the compiler happy.

}


// Convert a base to an array offset.
kgl::ContigOffset_t kgl::DNA5::nucleotideToColumn(Alphabet nucleotide) {

  // Translate the nucleotide to an array column
  switch (nucleotide) {

    case Alphabet::A: return A_NUCLEOTIDE_OFFSET;
    case Alphabet::C: return C_NUCLEOTIDE_OFFSET;
    case Alphabet::G: return G_NUCLEOTIDE_OFFSET;
    case Alphabet::T: return T_NUCLEOTIDE_OFFSET;
    case Alphabet::N: return N_NUCLEOTIDE_OFFSET;

  }

  return N_NUCLEOTIDE_OFFSET; //  Never reached, to keep the compiler happy.

}
