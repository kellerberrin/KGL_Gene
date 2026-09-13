//
// Created by kellerberrin on 17/11/17.
//


#include "kgl_alphabet_dna5.h"
#include "kgl_alphabet_dna5_tables.h"
#include "kel_exec_env.h"


namespace kgl = kellerberrin::genome;


namespace kellerberrin::genome::detail {   //  DNA5-specific tables (not shared with CodingDNA5).


// Extended IUPAC nucleotide codes (R,Y,S,W,K,M,B,D,H,V): accepted on input but converted to 'N'.
inline constexpr auto EXTENDED_TABLE = []() consteval {

  std::array<bool, 256> table{};
  table[static_cast<unsigned char>(DNA5::R_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::Y_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::S_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::W_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::K_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::M_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::B_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::D_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::H_NUCLEOTIDE)] = true;
  table[static_cast<unsigned char>(DNA5::V_NUCLEOTIDE)] = true;
  return table;

}();


// Complementary bases indexed over the unsigned char domain.
// Only the 5 canonical values map to a complement; every other (corrupted) value maps to N,
// matching the reference switch plus its trailing 'never reached' return.
inline constexpr auto COMPLEMENT_TABLE = []() consteval {

  std::array<CodingDNA5::Alphabet, 256> table{};
  table.fill(CodingDNA5::Alphabet::N);
  table[static_cast<unsigned char>(DNA5::Alphabet::A)] = CodingDNA5::Alphabet::T;
  table[static_cast<unsigned char>(DNA5::Alphabet::C)] = CodingDNA5::Alphabet::G;
  table[static_cast<unsigned char>(DNA5::Alphabet::G)] = CodingDNA5::Alphabet::C;
  table[static_cast<unsigned char>(DNA5::Alphabet::T)] = CodingDNA5::Alphabet::A;
  table[static_cast<unsigned char>(DNA5::Alphabet::N)] = CodingDNA5::Alphabet::N;
  return table;

}();


}   // end namespace


[[nodiscard]] bool kgl::DNA5::isExtended(Nucleotide_t char_letter) {

  return detail::EXTENDED_TABLE[static_cast<unsigned char>(char_letter)];

}


bool kgl::DNA5::validAlphabet(Alphabet nucleotide) {

  // We DO NOT use a switch here.
  // Because the switch assumes we can only have 5 base types (and a memory corrupted sequence may not).
  return detail::NucleotideTables<DNA5>::VALID_TABLE[static_cast<unsigned char>(nucleotide)];

}


// Convert char to Alphabet enum type.
kgl::DNA5::Alphabet kgl::DNA5::convertChar(char chr_base) {

  const Alphabet nucleotide = detail::NucleotideTables<DNA5>::ALPHABET_TABLE[static_cast<unsigned char>(chr_base)];

  // The reference switch distinguishes valid input characters (including 'N' in both cases,
  // which map to N) from the default error branch - the valid-character table preserves
  // exactly that distinction (A1).
  if (not detail::NucleotideTables<DNA5>::VALID_CHAR_TABLE[static_cast<unsigned char>(chr_base)]) {

    static bool report_extended = false;
    if (isExtended(chr_base) and not report_extended) {

      ExecEnv::log().warn("DNA5::convertChar(), IUPAC extended nucleotides detected, all converted to the unknown nucleotide 'N'");
      report_extended = true;

    } else if (not isExtended(chr_base)) {

      ExecEnv::log().error("DNA5::convertChar(), Unknown nucleotide detected: '{}', ascii value: {}. Input is probably corrupt or not DNA text.", chr_base, static_cast<size_t>(chr_base));

    }

  }

  return nucleotide;

}


// Find complementary bases.
kgl::CodingDNA5::Alphabet kgl::DNA5::complementNucleotide(Alphabet nucleotide) {

  // Total lookup over the unsigned char domain: a corrupted raw value cannot index
  // out of bounds and maps to N exactly as the reference did (A2).
  return detail::COMPLEMENT_TABLE[static_cast<unsigned char>(nucleotide)];

}


// Convert a base to an array offset.
kgl::ContigOffset_t kgl::DNA5::symbolToColumn(Alphabet nucleotide) {

  // The table lookup is total so a corrupted raw value cannot index out of bounds.
  return detail::NucleotideTables<DNA5>::COLUMN_TABLE[static_cast<unsigned char>(nucleotide)];

}


// Transitions involve interchanges of nucleotides of similar shapes: two-ring purines (A<>G)
// or one-ring pyrimidines (C<>T). A transversion_ is simply the complement of this function.
bool kgl::DNA5::isTransition(Alphabet nucleotide_1, Alphabet nucleotide_2) {

  return (nucleotide_1 == Alphabet::A and nucleotide_2 == Alphabet::G)
         or (nucleotide_1 == Alphabet::G and nucleotide_2 == Alphabet::A)
         or (nucleotide_1 == Alphabet::C and nucleotide_2 == Alphabet::T)
         or (nucleotide_1 == Alphabet::T and nucleotide_2 == Alphabet::C);

}


const std::vector<kgl::DNA5::Alphabet>& kgl::DNA5::enumerateAlphabet() {

  return detail::NucleotideTables<DNA5>::ALPHABET_VECTOR;

}