//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_ALPHABET_BASE_H
#define KGL_ALPHABET_BASE_H


#include <cstdint>
#include <memory>
#include <string>
#include <queue>
#include "kgl_logging.h"
#include "kgl_genome_types.h"
#include "kgl_exec_env.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


// Standard contig nucleotide column layout.
// Implement the NUCLEOTIDE_COLUMNS as "A", "C", "G", "T"/"U", "-", "E", "F", "I", "J" "K" in that order.
// Note that "E" = +A , "F" = +C , "I" = +G, "J" = +T/U and "K" = +N.


class DNA5 {

public:


  DNA5() = delete; // Singleton
  ~DNA5() = delete;

  static constexpr Nucleotide_t A_NUCLEOTIDE = 'A';
  static constexpr Nucleotide_t A_NUCLEOTIDE_LC = 'a';
  static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  static constexpr Nucleotide_t C_NUCLEOTIDE = 'C';
  static constexpr Nucleotide_t C_NUCLEOTIDE_LC = 'c';
  static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  static constexpr Nucleotide_t G_NUCLEOTIDE = 'G';
  static constexpr Nucleotide_t G_NUCLEOTIDE_LC = 'g';
  static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  static constexpr Nucleotide_t U_NUCLEOTIDE = 'U';
  static constexpr Nucleotide_t U_NUCLEOTIDE_LC = 'u';
  static constexpr ContigOffset_t U_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t T_NUCLEOTIDE = 'T';
  static constexpr Nucleotide_t T_NUCLEOTIDE_LC = 't';
  static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t N_NUCLEOTIDE = 'N';
  static constexpr Nucleotide_t N_NUCLEOTIDE_LC = 'n';
  static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;


  // The Alphabet enum type must be defined -see kgl_alphabet_string.h
  enum class Alphabet : Nucleotide_t
  { A = A_NUCLEOTIDE,
    C = C_NUCLEOTIDE,
    G = G_NUCLEOTIDE,
    T = T_NUCLEOTIDE,
    N = N_NUCLEOTIDE };

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  static Alphabet convertChar(char chr_base);

  // Find complementary bases.
  static Alphabet complementNucleotide(Alphabet nucleotide) {

    // Translate the nucleotide
    switch (nucleotide) {

      case Alphabet::A: return Alphabet::T;
      case Alphabet::C: return Alphabet::G;
      case Alphabet::G: return Alphabet::C;
      case Alphabet::T: return Alphabet::A;
      case Alphabet::N: return Alphabet::N;

    }

    return Alphabet::N; //  Never reached, to keep the compiler happy.

  }

  // Convert a base to an array offset.
  static ContigOffset_t nucleotideToColumn(Alphabet nucleotide) {

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

  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

};



class ExtendDNA5  {

public:


  ExtendDNA5() = delete; // Singleton

  static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 11;

  static constexpr Nucleotide_t A_NUCLEOTIDE = 'A';
  static constexpr Nucleotide_t A_NUCLEOTIDE_LC = 'a';
  static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  static constexpr Nucleotide_t C_NUCLEOTIDE = 'C';
  static constexpr Nucleotide_t C_NUCLEOTIDE_LC = 'c';
  static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  static constexpr Nucleotide_t G_NUCLEOTIDE = 'G';
  static constexpr Nucleotide_t G_NUCLEOTIDE_LC = 'g';
  static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  static constexpr Nucleotide_t U_NUCLEOTIDE = 'U';
  static constexpr Nucleotide_t U_NUCLEOTIDE_LC = 'u';
  static constexpr ContigOffset_t U_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t T_NUCLEOTIDE = 'T';
  static constexpr Nucleotide_t T_NUCLEOTIDE_LC = 't';
  static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t N_NUCLEOTIDE = 'N';
  static constexpr Nucleotide_t N_NUCLEOTIDE_LC = 'n';
  static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;
  static constexpr Nucleotide_t DELETE_NUCLEOTIDE = '-';
  static constexpr ContigOffset_t DELETE_NUCLEOTIDE_OFFSET = 5;
  static constexpr Nucleotide_t INSERT_A_NUCLEOTIDE = 'E';
  static constexpr ContigOffset_t INSERT_A_NUCLEOTIDE_OFFSET = 6;
  static constexpr Nucleotide_t INSERT_C_NUCLEOTIDE = 'F';
  static constexpr ContigOffset_t INSERT_C_NUCLEOTIDE_OFFSET = 7;
  static constexpr Nucleotide_t INSERT_G_NUCLEOTIDE = 'I';
  static constexpr ContigOffset_t INSERT_G_NUCLEOTIDE_OFFSET = 8;
  static constexpr Nucleotide_t INSERT_U_NUCLEOTIDE = 'J';
  static constexpr ContigOffset_t INSERT_U_NUCLEOTIDE_OFFSET = 9;
  static constexpr Nucleotide_t INSERT_T_NUCLEOTIDE = 'J';
  static constexpr ContigOffset_t INSERT_T_NUCLEOTIDE_OFFSET = 9;
  static constexpr Nucleotide_t INSERT_N_NUCLEOTIDE = 'K';
  static constexpr ContigOffset_t INSERT_N_NUCLEOTIDE_OFFSET = 10;


  // Note that X = Delete, E = A insert, F = C insert, I = G insert, J = T insert and K = N insert.
  // The Alphabet enum type must be defined -see kgl_alphabet_string.h
  // The extended alphabet does not conflict with LUPAC codes.
  enum class Alphabet : Nucleotide_t
  { A = A_NUCLEOTIDE,
    C = C_NUCLEOTIDE,
    G = G_NUCLEOTIDE,
    T = T_NUCLEOTIDE,
    N = N_NUCLEOTIDE,
    X = DELETE_NUCLEOTIDE,
    E = INSERT_A_NUCLEOTIDE,
    F = INSERT_C_NUCLEOTIDE,
    I = INSERT_G_NUCLEOTIDE,
    J = INSERT_T_NUCLEOTIDE,
    K = INSERT_N_NUCLEOTIDE
  };

  static bool isDeletion(Alphabet nucleotide) {

    return nucleotide == Alphabet::X;

  }

  static bool isInsertion(Alphabet nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case Alphabet::E:
      case Alphabet::F:
      case Alphabet::I:
      case Alphabet::J:
      case Alphabet::K:
        return true;

      default:
        return false;

    }

  }


  static bool isBaseCode(Alphabet nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case Alphabet::A:
      case Alphabet::C:
      case Alphabet::G:
      case Alphabet::T:
      case Alphabet::N:
        return true;

      default:
        return false;

    }

  }

  static DNA5::Alphabet extendToBase(Alphabet nucleotide) {

    switch (nucleotide) {

      case Alphabet::A: return DNA5::Alphabet::A;
      case Alphabet::C: return DNA5::Alphabet::C;
      case Alphabet::G: return DNA5::Alphabet::G;
      case Alphabet::T: return DNA5::Alphabet::T;
      case Alphabet::N: return DNA5::Alphabet::N;

      default:
        ExecEnv::log().error("ExtendDNA5::extendToBase, attempt to convert non-base nucleotide",
                             convertToChar(nucleotide));

    }

    return DNA5::Alphabet::N;

  }


  // Convert a base to an array offset.
  static ContigOffset_t nucleotideToColumn(Alphabet nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case Alphabet::A: return A_NUCLEOTIDE_OFFSET;
      case Alphabet::C: return C_NUCLEOTIDE_OFFSET;
      case Alphabet::G: return G_NUCLEOTIDE_OFFSET;
      case Alphabet::T: return T_NUCLEOTIDE_OFFSET;
      case Alphabet::N: return N_NUCLEOTIDE_OFFSET;
      case Alphabet::X: return DELETE_NUCLEOTIDE_OFFSET;
      case Alphabet::E: return INSERT_A_NUCLEOTIDE_OFFSET;
      case Alphabet::F: return INSERT_C_NUCLEOTIDE_OFFSET;
      case Alphabet::I: return INSERT_G_NUCLEOTIDE_OFFSET;
      case Alphabet::J: return INSERT_T_NUCLEOTIDE_OFFSET;
      case Alphabet::K: return INSERT_N_NUCLEOTIDE_OFFSET;

    }

    return N_NUCLEOTIDE_OFFSET; // Never reached, to keep the compiler happy.

  }

  // Converts an array offset into a base.
  static Alphabet offsetToNucleotide(ContigOffset_t offset) {

    // Translate the nucleotide to an array column
    switch (offset) {

      case A_NUCLEOTIDE_OFFSET: return Alphabet::A;
      case C_NUCLEOTIDE_OFFSET: return Alphabet::C;
      case G_NUCLEOTIDE_OFFSET: return Alphabet::G;
      case T_NUCLEOTIDE_OFFSET: return Alphabet::T;
      case N_NUCLEOTIDE_OFFSET: return Alphabet::N;
      case DELETE_NUCLEOTIDE_OFFSET: return Alphabet::X;
      case INSERT_A_NUCLEOTIDE_OFFSET: return Alphabet::E;
      case INSERT_C_NUCLEOTIDE_OFFSET: return Alphabet::F;
      case INSERT_G_NUCLEOTIDE_OFFSET: return Alphabet::I;
      case INSERT_T_NUCLEOTIDE_OFFSET: return Alphabet::J;
      case INSERT_N_NUCLEOTIDE_OFFSET: return Alphabet::K;

      default:
        ExecEnv::log().error("ExtendDNA5::offsetToNucleotide(), Invalid Nucleotide Offset", offset);
        return Alphabet::N;
    }

    return Alphabet::N; // Never reached, to keep the compiler happy.

  }


  // Insert offsets. Translates to an insert nucleotide.
  static Alphabet insertNucleotide(DNA5::Alphabet nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case DNA5::Alphabet::A: return Alphabet::E;
      case DNA5::Alphabet::C: return Alphabet::F;
      case DNA5::Alphabet::G: return Alphabet::I;
      case DNA5::Alphabet::T: return Alphabet::J;
      case DNA5::Alphabet::N: return Alphabet::K;

    }

    return Alphabet::K; // Dummy value

  }


  // Find extended complementary bases.
  static Alphabet complementNucleotide(Alphabet nucleotide) {

    switch (nucleotide) {

      case Alphabet::A: return Alphabet::T;
      case Alphabet::C: return Alphabet::G;
      case Alphabet::G: return Alphabet::C;
      case Alphabet::T: return Alphabet::A;
      case Alphabet::N: return Alphabet::N;
      case Alphabet::X: return Alphabet::X;
      case Alphabet::E: return Alphabet::J;
      case Alphabet::F: return Alphabet::I;
      case Alphabet::I: return Alphabet::F;
      case Alphabet::J: return Alphabet::E;
      case Alphabet::K: return Alphabet::K;

    }

    return Alphabet::N; // Never reached, to keep the compiler happy.

  }

  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  static Alphabet convertChar(char chr_base);

  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_ALPHABET_BASE_H
