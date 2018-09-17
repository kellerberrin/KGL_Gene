//
// Created by kellerberrin on 25/11/17.
//


#include "kgl_alphabet_coding_dna5.h"


namespace kgl = kellerberrin::genome;


const kgl::ContigOffset_t kgl::CodingDNA5::NUCLEOTIDE_COLUMNS;

const kgl::Nucleotide_t kgl::CodingDNA5::A_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::CodingDNA5::A_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::CodingDNA5::C_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::CodingDNA5::C_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::CodingDNA5::G_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::CodingDNA5::G_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::CodingDNA5::U_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::CodingDNA5::U_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::CodingDNA5::T_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::CodingDNA5::T_NUCLEOTIDE_OFFSET;
const kgl::Nucleotide_t kgl::CodingDNA5::N_NUCLEOTIDE;
const kgl::ContigOffset_t kgl::CodingDNA5::N_NUCLEOTIDE_OFFSET;


bool kgl::CodingDNA5::validAlphabet(Alphabet nucleotide) {

  auto int_value = static_cast<size_t>(nucleotide);

// We DO NOT use a switch here.
// Because the switch always assumes we can only have 5 base types (and a memory corrupted sequence may not).
  bool compare = int_value == static_cast<size_t>(Alphabet::A)
                 or int_value == static_cast<size_t>(Alphabet::C)
                 or int_value == static_cast<size_t>(Alphabet::G)
                 or int_value == static_cast<size_t>(Alphabet::T)
                 or int_value == static_cast<size_t>(Alphabet::N);

  return compare;

}


// Convert a base to an array offset.
kgl::ContigOffset_t kgl::CodingDNA5::nucleotideToColumn(Alphabet nucleotide) {

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
