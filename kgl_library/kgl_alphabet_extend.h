//
// Created by kellerberrin on 25/12/17.
//

#ifndef KGL_ALPHABET_EXTEND_H
#define KGL_ALPHABET_EXTEND_H


#include "kgl_logging.h"
#include "kgl_genome_types.h"
#include "kgl_exec_env.h"
#include "kgl_alphabet_dna5.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




// Implement the extended NUCLEOTIDE_COLUMNS as "A", "C", "G", "T"/"U", "-",  in that order.


class ExtendCountColumns  {

public:


  ExtendCountColumns() = delete; // Singleton

  static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 6;


  static constexpr Nucleotide_t A_NUCLEOTIDE = 'A';
  static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  static constexpr Nucleotide_t C_NUCLEOTIDE = 'C';
  static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  static constexpr Nucleotide_t G_NUCLEOTIDE = 'G';
  static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  static constexpr Nucleotide_t U_NUCLEOTIDE = 'U';
  static constexpr ContigOffset_t U_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t T_NUCLEOTIDE = 'T';
  static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_t N_NUCLEOTIDE = 'N';
  static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;
  static constexpr Nucleotide_t DELETE_NUCLEOTIDE = '-';
  static constexpr ContigOffset_t DELETE_NUCLEOTIDE_OFFSET = 5;


  // Note that X = Delete
  // The extended alphabet does not conflict with LUPAC codes.
  enum class Alphabet : Nucleotide_t
  { A = A_NUCLEOTIDE,
    C = C_NUCLEOTIDE,
    G = G_NUCLEOTIDE,
    T = T_NUCLEOTIDE,
    N = N_NUCLEOTIDE,
    X = DELETE_NUCLEOTIDE,
  };

  static bool isDeletion(Alphabet nucleotide) {

    return nucleotide == Alphabet::X;

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
        ExecEnv::log().error("ExtendCountColumns::extendToBase, attempt to convert non-base nucleotide",
                             convertToChar(nucleotide));

    }

    return DNA5::Alphabet::N;

  }

  static Alphabet baseToExtend(DNA5::Alphabet nucleotide) {

    switch (nucleotide) {

      case DNA5::Alphabet::A: return Alphabet::A;
      case DNA5::Alphabet::C: return Alphabet::C;
      case DNA5::Alphabet::G: return Alphabet::G;
      case DNA5::Alphabet::T: return Alphabet::T;
      case DNA5::Alphabet::N: return Alphabet::N;

    }

    return Alphabet::N; // Never reached.

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

      default:
        ExecEnv::log().error("ReadCountColumns::offsetToNucleotide(), Invalid Nucleotide Offset", offset);
        return Alphabet::N;
    }

    return Alphabet::N; // Never reached, to keep the compiler happy.

  }


  // The Alphabet convertChar(char) function must be defined -see kgl_alphabet_string.h
  static Alphabet convertChar(char chr_base);

  static char convertToChar(Alphabet nucleotide) { return static_cast<char>(nucleotide); }

};



}   // namespace genome
}   // namespace kellerberrin





#endif //KGL_ALPHABET_EXTEND_H
