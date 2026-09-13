//
// Created by kellerberrin on 17/11/17.
//

#ifndef KGL_ALPHABET_DNA5_TABLES_H
#define KGL_ALPHABET_DNA5_TABLES_H


// Internal header: shared constexpr lookup tables for the DNA5 and CodingDNA5 alphabets.
// Not part of the public interface; include only from the alphabet translation units.

#include <array>
#include <vector>

#include "kgl_alphabet_coding_dna5.h"
#include "kgl_alphabet_dna5.h"


namespace kellerberrin::genome::detail {   //  Shared implementation for the DNA5 alphabets (single-sourced).


// The constexpr lookup tables are total over the unsigned char domain so the
// conversion functions are branch-free and a memory-corrupted raw value can
// never index out of bounds (the no-switch intent of the original code).


template<typename AlphabetType>
struct NucleotideTables {

  /// Valid input characters (both cases; A/C/G/T/U/N) - used by the convertChar error branch.
  inline static constexpr auto VALID_CHAR_TABLE = []() consteval {

    std::array<bool, 256> table{};
    for (const bool lower_case : {false, true}) {

      const char case_offset = lower_case ? ('a' - 'A') : 0;
      table[static_cast<unsigned char>(AlphabetType::A_NUCLEOTIDE + case_offset)] = true;
      table[static_cast<unsigned char>(AlphabetType::C_NUCLEOTIDE + case_offset)] = true;
      table[static_cast<unsigned char>(AlphabetType::G_NUCLEOTIDE + case_offset)] = true;
      table[static_cast<unsigned char>(AlphabetType::U_NUCLEOTIDE + case_offset)] = true;
      table[static_cast<unsigned char>(AlphabetType::T_NUCLEOTIDE + case_offset)] = true;
      table[static_cast<unsigned char>(AlphabetType::N_NUCLEOTIDE + case_offset)] = true;

    }
    return table;

  }();

  /// Valid alphabet enum values. The enum defines only A/C/G/T/N; a raw 'U' byte in an
  /// enum slot is a corrupted value and is invalid - exactly as the reference or-chain was.
  inline static constexpr auto VALID_TABLE = []() consteval {

    std::array<bool, 256> table{};
    table[static_cast<unsigned char>(AlphabetType::A_NUCLEOTIDE)] = true;
    table[static_cast<unsigned char>(AlphabetType::C_NUCLEOTIDE)] = true;
    table[static_cast<unsigned char>(AlphabetType::G_NUCLEOTIDE)] = true;
    table[static_cast<unsigned char>(AlphabetType::T_NUCLEOTIDE)] = true;
    table[static_cast<unsigned char>(AlphabetType::N_NUCLEOTIDE)] = true;
    return table;

  }();

  /// Alphabet enum value for each input character (both cases; U maps to T as the reference did).
  inline static constexpr auto ALPHABET_TABLE = []() consteval {

    using Alphabet = typename AlphabetType::Alphabet;
    std::array<Alphabet, 256> table{};
    table.fill(Alphabet::N);
    // Both cases are mapped (the reference switch upper-cased via std::toupper).
    for (const bool lower_case : {false, true}) {

      const char case_offset = lower_case ? ('a' - 'A') : 0;
      table[static_cast<unsigned char>(AlphabetType::A_NUCLEOTIDE + case_offset)] = Alphabet::A;
      table[static_cast<unsigned char>(AlphabetType::C_NUCLEOTIDE + case_offset)] = Alphabet::C;
      table[static_cast<unsigned char>(AlphabetType::G_NUCLEOTIDE + case_offset)] = Alphabet::G;
      table[static_cast<unsigned char>(AlphabetType::U_NUCLEOTIDE + case_offset)] = Alphabet::T;
      table[static_cast<unsigned char>(AlphabetType::T_NUCLEOTIDE + case_offset)] = Alphabet::T;
      table[static_cast<unsigned char>(AlphabetType::N_NUCLEOTIDE + case_offset)] = Alphabet::N;

    }
    return table;

  }();

  /// Column offset for each alphabet enum value. The reference switch enumerated only
  /// A/C/G/T/N; a raw 'U' byte in an enum slot is corrupted and maps to the N column.
  inline static constexpr auto COLUMN_TABLE = []() consteval {

    std::array<ContigOffset_t, 256> table{};
    table.fill(AlphabetType::N_NUCLEOTIDE_OFFSET);
    table[static_cast<unsigned char>(AlphabetType::A_NUCLEOTIDE)] = AlphabetType::A_NUCLEOTIDE_OFFSET;
    table[static_cast<unsigned char>(AlphabetType::C_NUCLEOTIDE)] = AlphabetType::C_NUCLEOTIDE_OFFSET;
    table[static_cast<unsigned char>(AlphabetType::G_NUCLEOTIDE)] = AlphabetType::G_NUCLEOTIDE_OFFSET;
    table[static_cast<unsigned char>(AlphabetType::T_NUCLEOTIDE)] = AlphabetType::T_NUCLEOTIDE_OFFSET;
    table[static_cast<unsigned char>(AlphabetType::N_NUCLEOTIDE)] = AlphabetType::N_NUCLEOTIDE_OFFSET;
    return table;

  }();

  /// The enumerated alphabet vector (A, C, G, T, N).
  inline static const std::vector<typename AlphabetType::Alphabet>& ALPHABET_VECTOR = []() {

    using Alphabet = typename AlphabetType::Alphabet;
    static std::vector<Alphabet> alphabet_vector = {Alphabet::A, Alphabet::C, Alphabet::G, Alphabet::T, Alphabet::N};
    return alphabet_vector;

  }();

};


}   // end namespace


#endif //KGL_ALPHABET_DNA5_TABLES_H