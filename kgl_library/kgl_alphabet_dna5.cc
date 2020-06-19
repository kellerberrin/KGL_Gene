//
// Created by kellerberrin on 17/11/17.
//


#include "kgl_alphabet_dna5.h"


namespace kgl = kellerberrin::genome;



[[nodiscard]] bool kgl::DNA5::isExtended(Nucleotide_t char_letter) {

  switch(char_letter) {

    case R_NUCLEOTIDE:
    case Y_NUCLEOTIDE:
    case S_NUCLEOTIDE:
    case W_NUCLEOTIDE:
    case K_NUCLEOTIDE:
    case M_NUCLEOTIDE:
    case B_NUCLEOTIDE:
    case D_NUCLEOTIDE:
    case H_NUCLEOTIDE:
    case V_NUCLEOTIDE:
      return true;

    default:
      return false;
  }

}


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

    default: {

      static bool report_extended = false;
      if (isExtended(chr_base) && not report_extended) {

        ExecEnv::log().warn("DNA5::convertchar(), IUPAC extended nucleotides detected, all converted to the unknown nucleotide 'N'");
        report_extended = true;

      } else if (not isExtended(chr_base)) {

        ExecEnv::log().error("DNA5::convertchar(), Unknown nucleotide detected: '{}', ascii value: {}. Input is probably corrupt or not DNA text.", chr_base, static_cast<size_t>(chr_base));

      }
      return Alphabet::N;

    }

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
kgl::ContigOffset_t kgl::DNA5::symbolToColumn(Alphabet nucleotide) {

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

// Transitions involve interchanges of nucleotides of similar shapes: two-ring purines (A<>G)
// or one-ring pyrimidines (C<>T). A transversion is simply the complement of this function.
bool kgl::DNA5::isTransition(Alphabet nucleotide_1, Alphabet nucleotide_2) {

  return (nucleotide_1 == Alphabet::A and nucleotide_2 == Alphabet::G)
         or (nucleotide_1 == Alphabet::G and nucleotide_2 == Alphabet::A)
         or (nucleotide_1 == Alphabet::C and nucleotide_2 == Alphabet::T)
         or (nucleotide_1 == Alphabet::T and nucleotide_2 == Alphabet::C);

}


const std::vector<kgl::DNA5::Alphabet>& kgl::DNA5::enumerateAlphabet() {

  static std::vector<Alphabet> alphabet_vector = {Alphabet::A, Alphabet::C, Alphabet::G, Alphabet::T, Alphabet::N};

  return alphabet_vector;

}
