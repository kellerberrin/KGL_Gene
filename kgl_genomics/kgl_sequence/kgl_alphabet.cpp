//
// kgl_alphabet.cpp — corruption-safe total lookup tables for the alphabet policies.
//
// The tables are total over the unsigned char domain, so a memory-corrupted raw value can
// never index out of bounds (the no-switch intent of the original code).
//

#include "kgl_alphabet.h"
#include "kel_exec_env.h"

#include <array>
#include <vector>


namespace kgl = kellerberrin::genome;


namespace kellerberrin::genome::detail {   //  Single-sourced nucleotide tables (DNA5 + CodingDNA5).


// Valid input characters (both cases; A/C/G/T/U/N) — used by the convertChar error branch.
inline constexpr auto VALID_CHAR_TABLE = []() consteval {

  std::array<bool, 256> table{};
  for (const bool lower_case : {false, true}) {

    const char case_offset = lower_case ? ('a' - 'A') : 0;
    table[static_cast<unsigned char>('A' + case_offset)] = true;
    table[static_cast<unsigned char>('C' + case_offset)] = true;
    table[static_cast<unsigned char>('G' + case_offset)] = true;
    table[static_cast<unsigned char>('U' + case_offset)] = true;
    table[static_cast<unsigned char>('T' + case_offset)] = true;
    table[static_cast<unsigned char>('N' + case_offset)] = true;

  }
  return table;

}();


// Valid alphabet enum values. The enum defines only A/C/G/T/N; a raw 'U' byte in an
// enum slot is a corrupted value and is invalid — exactly as the reference or-chain was.
inline constexpr auto VALID_TABLE = []() consteval {

  std::array<bool, 256> table{};
  table[static_cast<unsigned char>(Nucleotide::A)] = true;
  table[static_cast<unsigned char>(Nucleotide::C)] = true;
  table[static_cast<unsigned char>(Nucleotide::G)] = true;
  table[static_cast<unsigned char>(Nucleotide::T)] = true;
  table[static_cast<unsigned char>(Nucleotide::N)] = true;
  return table;

}();


// Alphabet enum value for each input character (both cases; U maps to T as the reference did).
inline constexpr auto ALPHABET_TABLE = []() consteval {

  std::array<Nucleotide, 256> table{};
  table.fill(Nucleotide::N);
  for (const bool lower_case : {false, true}) {

    const char case_offset = lower_case ? ('a' - 'A') : 0;
    table[static_cast<unsigned char>('A' + case_offset)] = Nucleotide::A;
    table[static_cast<unsigned char>('C' + case_offset)] = Nucleotide::C;
    table[static_cast<unsigned char>('G' + case_offset)] = Nucleotide::G;
    table[static_cast<unsigned char>('U' + case_offset)] = Nucleotide::T;
    table[static_cast<unsigned char>('T' + case_offset)] = Nucleotide::T;
    table[static_cast<unsigned char>('N' + case_offset)] = Nucleotide::N;

  }
  return table;

}();


// Column offset for each alphabet enum value. A raw 'U' enum byte is corrupted and maps to N.
inline constexpr auto COLUMN_TABLE = []() consteval {

  std::array<ContigOffset_t, 256> table{};
  table.fill(4);   // N column default (corruption-safe)
  table[static_cast<unsigned char>(Nucleotide::A)] = 0;
  table[static_cast<unsigned char>(Nucleotide::C)] = 1;
  table[static_cast<unsigned char>(Nucleotide::G)] = 2;
  table[static_cast<unsigned char>(Nucleotide::T)] = 3;
  table[static_cast<unsigned char>(Nucleotide::N)] = 4;
  return table;

}();


// Extended IUPAC nucleotide codes (R,Y,S,W,K,M,B,D,H,V): accepted on input but converted to 'N'.
inline constexpr auto EXTENDED_TABLE = []() consteval {

  std::array<bool, 256> table{};
  for (const char code : {'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}) {
    table[static_cast<unsigned char>(code)] = true;
  }
  return table;

}();


// Complementary bases. Only the 5 canonical values map to a complement; every corrupted
// value maps to N (matching the reference switch plus its trailing 'never reached' return).
inline constexpr auto COMPLEMENT_TABLE = []() consteval {

  std::array<Nucleotide, 256> table{};
  table.fill(Nucleotide::N);
  table[static_cast<unsigned char>(Nucleotide::A)] = Nucleotide::T;
  table[static_cast<unsigned char>(Nucleotide::C)] = Nucleotide::G;
  table[static_cast<unsigned char>(Nucleotide::G)] = Nucleotide::C;
  table[static_cast<unsigned char>(Nucleotide::T)] = Nucleotide::A;
  table[static_cast<unsigned char>(Nucleotide::N)] = Nucleotide::N;
  return table;

}();


// The enumerated nucleotide alphabet (A, C, G, T, N).
inline constexpr std::array<Nucleotide, 5> NUCLEOTIDE_ALPHABET{Nucleotide::A, Nucleotide::C,
                                                               Nucleotide::G, Nucleotide::T,
                                                               Nucleotide::N};


}   // end namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA5
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::DNA5::validAlphabet(Alphabet nucleotide) noexcept {

  return detail::VALID_TABLE[static_cast<unsigned char>(nucleotide)];

}


kgl::DNA5::Alphabet kgl::DNA5::convertChar(char chr_base) {

  // Silent and total: any invalid char maps to N. Diagnostics are reported once per parsed
  // string by Sequence's parse constructor via reportInvalidDNA5() (improvement over the
  // reference's per-character log with a racy function-local static bool).
  return detail::ALPHABET_TABLE[static_cast<unsigned char>(chr_base)];

}


bool kgl::DNA5::isExtended(char char_letter) noexcept {

  return detail::EXTENDED_TABLE[static_cast<unsigned char>(char_letter)];

}


kgl::Nucleotide kgl::DNA5::complementNucleotide(Alphabet nucleotide) noexcept {

  return detail::COMPLEMENT_TABLE[static_cast<unsigned char>(nucleotide)];

}


kgl::ContigOffset_t kgl::DNA5::symbolToColumn(Alphabet nucleotide) noexcept {

  return detail::COLUMN_TABLE[static_cast<unsigned char>(nucleotide)];

}


bool kgl::DNA5::isTransition(Alphabet nucleotide_1, Alphabet nucleotide_2) noexcept {

  return (nucleotide_1 == Alphabet::A and nucleotide_2 == Alphabet::G)
         or (nucleotide_1 == Alphabet::G and nucleotide_2 == Alphabet::A)
         or (nucleotide_1 == Alphabet::C and nucleotide_2 == Alphabet::T)
         or (nucleotide_1 == Alphabet::T and nucleotide_2 == Alphabet::C);

}


const std::array<kgl::DNA5::Alphabet, 5>& kgl::DNA5::enumerateAlphabet() noexcept {

  return detail::NUCLEOTIDE_ALPHABET;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CodingDNA5 — shares the nucleotide implementation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::CodingDNA5::validAlphabet(Alphabet nucleotide) noexcept {

  return detail::VALID_TABLE[static_cast<unsigned char>(nucleotide)];

}


kgl::ContigOffset_t kgl::CodingDNA5::symbolToColumn(Alphabet nucleotide) noexcept {

  return detail::COLUMN_TABLE[static_cast<unsigned char>(nucleotide)];

}


const std::array<kgl::CodingDNA5::Alphabet, 5>& kgl::CodingDNA5::enumerateAlphabet() noexcept {

  return detail::NUCLEOTIDE_ALPHABET;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AminoAcid
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin::genome::detail {   //  Amino acid tables.


inline constexpr auto AMINO_COLUMN_TABLE = []() consteval {

  std::array<ContigOffset_t, 256> table{};
  table.fill(21);   // unknown amino column
  table[static_cast<unsigned char>(Amino::F)] = 0;
  table[static_cast<unsigned char>(Amino::L)] = 1;
  table[static_cast<unsigned char>(Amino::S)] = 2;
  table[static_cast<unsigned char>(Amino::Y)] = 3;
  table[static_cast<unsigned char>(Amino::C)] = 4;
  table[static_cast<unsigned char>(Amino::W)] = 5;
  table[static_cast<unsigned char>(Amino::P)] = 6;
  table[static_cast<unsigned char>(Amino::H)] = 7;
  table[static_cast<unsigned char>(Amino::Q)] = 8;
  table[static_cast<unsigned char>(Amino::R)] = 9;
  table[static_cast<unsigned char>(Amino::I)] = 10;
  table[static_cast<unsigned char>(Amino::M)] = 11;
  table[static_cast<unsigned char>(Amino::T)] = 12;
  table[static_cast<unsigned char>(Amino::N)] = 13;
  table[static_cast<unsigned char>(Amino::K)] = 14;
  table[static_cast<unsigned char>(Amino::V)] = 15;
  table[static_cast<unsigned char>(Amino::A)] = 16;
  table[static_cast<unsigned char>(Amino::D)] = 17;
  table[static_cast<unsigned char>(Amino::E)] = 18;
  table[static_cast<unsigned char>(Amino::G)] = 19;
  table[static_cast<unsigned char>(Amino::Stop)] = 20;
  table[static_cast<unsigned char>(Amino::Z)] = 21;
  return table;

}();


// Includes the rare selenocysteine (U) and pyrrolysine (O) so a parsed sequence containing
// them passes the module's own verifyString() corruption check.
inline constexpr auto AMINO_VALID_TABLE = []() consteval {

  std::array<bool, 256> table{};
  for (std::size_t index = 0; index < 256; ++index) {

    table[index] = (index == static_cast<unsigned char>(Amino::Z))
                   or (AMINO_COLUMN_TABLE[index] != 21);

  }
  table[static_cast<unsigned char>(Amino::U)] = true;
  table[static_cast<unsigned char>(Amino::O)] = true;
  return table;

}();


inline constexpr auto AMINO_ALPHABET_TABLE = []() consteval {

  std::array<Amino, 256> table{};
  table.fill(Amino::Z);
  table[static_cast<unsigned char>(Amino::F)] = Amino::F;
  table[static_cast<unsigned char>(Amino::L)] = Amino::L;
  table[static_cast<unsigned char>(Amino::S)] = Amino::S;
  table[static_cast<unsigned char>(Amino::Y)] = Amino::Y;
  table[static_cast<unsigned char>(Amino::C)] = Amino::C;
  table[static_cast<unsigned char>(Amino::W)] = Amino::W;
  table[static_cast<unsigned char>(Amino::P)] = Amino::P;
  table[static_cast<unsigned char>(Amino::H)] = Amino::H;
  table[static_cast<unsigned char>(Amino::Q)] = Amino::Q;
  table[static_cast<unsigned char>(Amino::R)] = Amino::R;
  table[static_cast<unsigned char>(Amino::I)] = Amino::I;
  table[static_cast<unsigned char>(Amino::M)] = Amino::M;
  table[static_cast<unsigned char>(Amino::T)] = Amino::T;
  table[static_cast<unsigned char>(Amino::N)] = Amino::N;
  table[static_cast<unsigned char>(Amino::K)] = Amino::K;
  table[static_cast<unsigned char>(Amino::V)] = Amino::V;
  table[static_cast<unsigned char>(Amino::A)] = Amino::A;
  table[static_cast<unsigned char>(Amino::D)] = Amino::D;
  table[static_cast<unsigned char>(Amino::E)] = Amino::E;
  table[static_cast<unsigned char>(Amino::G)] = Amino::G;
  table[static_cast<unsigned char>(Amino::U)] = Amino::U;
  table[static_cast<unsigned char>(Amino::O)] = Amino::O;
  table[static_cast<unsigned char>(Amino::Stop)] = Amino::Stop;
  table[static_cast<unsigned char>(Amino::Z)] = Amino::Z;
  return table;

}();


inline constexpr std::array<Amino, 22> AMINO_ALPHABET{Amino::F, Amino::L, Amino::S, Amino::Y,
                                                       Amino::C, Amino::W, Amino::P, Amino::H,
                                                       Amino::Q, Amino::R, Amino::I, Amino::M,
                                                       Amino::T, Amino::N, Amino::K, Amino::V,
                                                       Amino::A, Amino::D, Amino::E, Amino::G,
                                                       Amino::Stop, Amino::Z};


}   // end namespace


kgl::AminoAcid::Alphabet kgl::AminoAcid::convertChar(char chr_aa) {

  // Silent and total (strict: lowercase is NOT accepted, matching the reference). Diagnostics
  // are reported once per parsed string by Sequence's parse constructor via
  // reportInvalidAminoAcid() (improvement over per-character logging).
  return detail::AMINO_ALPHABET_TABLE[static_cast<unsigned char>(chr_aa)];

}


kgl::ContigOffset_t kgl::AminoAcid::symbolToColumn(Alphabet amino) noexcept {

  const ContigOffset_t column = detail::AMINO_COLUMN_TABLE[static_cast<unsigned char>(amino)];

  // Restore the reference corruption diagnostic (valid U/O are silent).
  if (not validAlphabet(amino)) {

    ExecEnv::log().error("AminoAcid::symbolToColumn(), Invalid amino symbol: {}", static_cast<char>(amino));

  }

  return column;

}


bool kgl::AminoAcid::validAlphabet(Alphabet amino) noexcept {

  return detail::AMINO_VALID_TABLE[static_cast<unsigned char>(amino)];

}


const std::array<kgl::AminoAcid::Alphabet, 22>& kgl::AminoAcid::enumerateAlphabet() noexcept {

  return detail::AMINO_ALPHABET;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Aggregate parse diagnostics.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin::genome {

[[nodiscard]] ParseReport reportInvalidDNA5(std::string_view alphabet_str) noexcept {

  ParseReport report;
  for (const char chr_base : alphabet_str) {
    detail::tallyInvalidDNA5(report, chr_base);
  }
  return report;

}


[[nodiscard]] ParseReport reportInvalidAminoAcid(std::string_view alphabet_str) noexcept {

  ParseReport report;
  for (const char chr_aa : alphabet_str) {
    detail::tallyInvalidAminoAcid(report, chr_aa);
  }
  return report;

}

}   // end namespace


namespace kellerberrin::genome::detail {

void tallyInvalidDNA5(ParseReport& report, char chr_base) noexcept {

  const auto raw = static_cast<unsigned char>(chr_base);
  if (not VALID_CHAR_TABLE[raw]) {

    if (DNA5::isExtended(chr_base)) {
      ++report.extended_chars;
    } else {
      ++report.invalid_chars;
    }

  }

}


void tallyInvalidAminoAcid(ParseReport& report, char chr_aa) noexcept {

  if (AMINO_ALPHABET_TABLE[static_cast<unsigned char>(chr_aa)] == Amino::Z
      and chr_aa != AminoAcid::UNKNOWN_AMINO) {

    ++report.invalid_chars;

  }

}

}   // end namespace
