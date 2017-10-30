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
// Implement the NUCLEOTIDE_COLUMNS as "A", "C", "G", "T"/"U", "-", "+" in that order.

class NucleotideColumn_DNA5 {

public:

  using NucleotideType = Nucleotide_DNA5_t;

  NucleotideColumn_DNA5() = delete; // Singleton
  ~NucleotideColumn_DNA5() = delete;

  static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 7;
  static constexpr Nucleotide_DNA5_t A_NUCLEOTIDE = 'A';
  static constexpr Nucleotide_DNA5_t A_NUCLEOTIDE_LC = 'a';
  static constexpr ContigOffset_t A_NUCLEOTIDE_OFFSET = 0;
  static constexpr Nucleotide_DNA5_t C_NUCLEOTIDE = 'C';
  static constexpr Nucleotide_DNA5_t C_NUCLEOTIDE_LC = 'c';
  static constexpr ContigOffset_t C_NUCLEOTIDE_OFFSET = 1;
  static constexpr Nucleotide_DNA5_t G_NUCLEOTIDE = 'G';
  static constexpr Nucleotide_DNA5_t G_NUCLEOTIDE_LC = 'g';
  static constexpr ContigOffset_t G_NUCLEOTIDE_OFFSET = 2;
  static constexpr Nucleotide_DNA5_t U_NUCLEOTIDE = 'U';
  static constexpr Nucleotide_DNA5_t U_NUCLEOTIDE_LC = 'u';
  static constexpr ContigOffset_t U_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_DNA5_t T_NUCLEOTIDE = 'T';
  static constexpr Nucleotide_DNA5_t T_NUCLEOTIDE_LC = 't';
  static constexpr ContigOffset_t T_NUCLEOTIDE_OFFSET = 3;
  static constexpr Nucleotide_DNA5_t N_NUCLEOTIDE = 'N';
  static constexpr Nucleotide_DNA5_t N_NUCLEOTIDE_LC = 'n';
  static constexpr ContigOffset_t N_NUCLEOTIDE_OFFSET = 4;
  static constexpr Nucleotide_DNA5_t DELETE_NUCLEOTIDE = '-';
  static constexpr ContigOffset_t DELETE_NUCLEOTIDE_OFFSET = 5;
  static constexpr Nucleotide_DNA5_t INSERT_SEQUENCE = '+';
  static constexpr ContigOffset_t INSERT_NUCLEOTIDE_OFFSET = 6;

  // Convert a base to an array offset.
  static ContigOffset_t nucleotideToColumn(const Nucleotide_DNA5_t nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case A_NUCLEOTIDE:
      case A_NUCLEOTIDE_LC: return A_NUCLEOTIDE_OFFSET;

      case C_NUCLEOTIDE:
      case C_NUCLEOTIDE_LC: return C_NUCLEOTIDE_OFFSET;

      case G_NUCLEOTIDE:
      case G_NUCLEOTIDE_LC: return G_NUCLEOTIDE_OFFSET;

      case T_NUCLEOTIDE:
      case T_NUCLEOTIDE_LC:
      case U_NUCLEOTIDE:
      case U_NUCLEOTIDE_LC: return T_NUCLEOTIDE_OFFSET;

      case N_NUCLEOTIDE:
      case N_NUCLEOTIDE_LC: return N_NUCLEOTIDE_OFFSET;

      case DELETE_NUCLEOTIDE: return DELETE_NUCLEOTIDE_OFFSET;

      case INSERT_SEQUENCE: return INSERT_NUCLEOTIDE_OFFSET;

      default:
        ExecEnv::log().critical("nucleotideToColumn(), Nucleotide array accessed with unknown nucleotide: {}",
                                nucleotide);
        return 0; // Never reached, to keep the compiler happy.

    }

  }

  // Converts an array offset into a base.
  static Nucleotide_DNA5_t offsetToNucleotide(ContigOffset_t offset) {

    // Translate the nucleotide to an array column
    switch (offset) {

      case A_NUCLEOTIDE_OFFSET:
        return A_NUCLEOTIDE;

      case C_NUCLEOTIDE_OFFSET:
        return C_NUCLEOTIDE;

      case G_NUCLEOTIDE_OFFSET:
        return G_NUCLEOTIDE;

      case T_NUCLEOTIDE_OFFSET:
        return T_NUCLEOTIDE;

      case N_NUCLEOTIDE_OFFSET:
        return N_NUCLEOTIDE;

      case DELETE_NUCLEOTIDE_OFFSET:
        return DELETE_NUCLEOTIDE;

      case INSERT_NUCLEOTIDE_OFFSET:
        return INSERT_SEQUENCE;

      default:
        ExecEnv::log().critical("indexToNucleotide(), unknown nucleotide offset: {}", offset);
        return N_NUCLEOTIDE; // Never reached, to keep the compiler happy.

    }

  }

  // Find complementary bases.
  static Nucleotide_DNA5_t complementNucleotide(const Nucleotide_DNA5_t nucleotide) {

    // Translate the nucleotide to an array column
    switch (nucleotide) {

      case A_NUCLEOTIDE:
      case A_NUCLEOTIDE_LC:
        return T_NUCLEOTIDE;

      case C_NUCLEOTIDE:
      case C_NUCLEOTIDE_LC:
        return G_NUCLEOTIDE;

      case G_NUCLEOTIDE:
      case G_NUCLEOTIDE_LC:
        return C_NUCLEOTIDE;

      case T_NUCLEOTIDE:
      case T_NUCLEOTIDE_LC:
        return A_NUCLEOTIDE;

      case U_NUCLEOTIDE:
      case U_NUCLEOTIDE_LC:
        return A_NUCLEOTIDE;

      default:
        ExecEnv::log().error("complementNucleotide(), Unknown Nucleotide: {} (ascii): {}",
                             nucleotide, static_cast<unsigned char>(nucleotide));
        return N_NUCLEOTIDE; // Never reached, to keep the compiler happy.

    }

  }

private:

};

}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_ALPHABET_BASE_H
