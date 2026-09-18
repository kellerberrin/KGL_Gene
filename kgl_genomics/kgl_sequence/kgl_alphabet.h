//
// kgl_alphabet.h — refactored alphabet policies (deepseek-v4-refactor).
//
// Single symbol domain per chemistry (Nucleotide, Amino) plus policy structs (DNA5,
// CodingDNA5, AminoAcid) that retain the reference's static "Alphabet::" protocol.
// The corruption-safe total lookup tables live in the .cpp translation unit.
//

#ifndef KGL_ALPHABET_H
#define KGL_ALPHABET_H


#include <array>
#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

#include "kgl_genome_types.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The symbol domains.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class Nucleotide : char { A = 'A', C = 'C', G = 'G', T = 'T', N = 'N' };

enum class Amino : char { F = 'F', L = 'L', S = 'S', Y = 'Y', C = 'C', W = 'W', P = 'P', H = 'H',
                          Q = 'Q', R = 'R', I = 'I', M = 'M', T = 'T', N = 'N', K = 'K', V = 'V',
                          A = 'A', D = 'D', E = 'E', G = 'G', U = 'U', O = 'O',
                          Stop = '*', Z = 'Z' };

static_assert(sizeof(Nucleotide) == 1 and sizeof(Amino) == 1, "Alphabet symbols must be byte sized");


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ParseReport — aggregate parse diagnostics (improvement over the reference's per-character
// logging with a racy function-local static bool). Counts; does not log. Sequence's parse
// constructor reports one aggregate diagnostic per parsed string.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ParseReport {

  std::size_t invalid_chars{0};    // Chars not in the alphabet and not extended IUPAC.
  std::size_t extended_chars{0};   // Extended IUPAC codes (DNA only; all converted to N).

};

/// Scan an ASCII string, counting chars that map to the unknown symbol. Thread-safe.
[[nodiscard]] ParseReport reportInvalidDNA5(std::string_view alphabet_str) noexcept;
/// Scan a protein string, counting chars that map to the unknown amino. Thread-safe.
[[nodiscard]] ParseReport reportInvalidAminoAcid(std::string_view alphabet_str) noexcept;

namespace detail {

  // Single-character tallies over the total lookup tables. Sequence's parse constructor calls
  // these inside its conversion loop, so a parsed string is scanned ONCE (v4 improvement over
  // the separate report pass). The public reportInvalid* wrappers loop over these helpers.
  void tallyInvalidDNA5(ParseReport& report, char chr_base) noexcept;
  void tallyInvalidAminoAcid(ParseReport& report, char chr_aa) noexcept;

}   // namespace detail


// constexpr char -> amino symbol (no logging; safe to use in consteval table decoding).
[[nodiscard]] constexpr Amino aminoFromChar(char amino_char) noexcept {

  switch (amino_char) {
    case 'F': return Amino::F;  case 'L': return Amino::L;  case 'S': return Amino::S;
    case 'Y': return Amino::Y;  case 'C': return Amino::C;  case 'W': return Amino::W;
    case 'P': return Amino::P;  case 'H': return Amino::H;  case 'Q': return Amino::Q;
    case 'R': return Amino::R;  case 'I': return Amino::I;  case 'M': return Amino::M;
    case 'T': return Amino::T;  case 'N': return Amino::N;  case 'K': return Amino::K;
    case 'V': return Amino::V;  case 'A': return Amino::A;  case 'D': return Amino::D;
    case 'E': return Amino::E;  case 'G': return Amino::G;  case 'U': return Amino::U;
    case 'O': return Amino::O;  case '*': return Amino::Stop;
    default:  return Amino::Z;
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA5 — the unstranded DNA policy. Retains the reference's static interface names.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct DNA5 {

  using Alphabet = Nucleotide;

  inline static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 5;

  DNA5() = delete;   // policy namespace holder
  ~DNA5() = delete;

  /// ReturnType a vector of all valid alphabet values.
  [[nodiscard]] static const std::array<Alphabet, 5>& enumerateAlphabet() noexcept;

  /// Checks for possible memory corruption.
  [[nodiscard]] static bool validAlphabet(Alphabet nucleotide) noexcept;

  /// Convert char to Alphabet enum type (logs on invalid / extended input).
  [[nodiscard]] static Alphabet convertChar(char chr_base);

  /// Extended IUPAC codes are recognised but converted to 'N' on input.
  [[nodiscard]] static bool isExtended(char char_letter) noexcept;

  /// Find complementary bases.
  [[nodiscard]] static Nucleotide complementNucleotide(Alphabet nucleotide) noexcept;

  /// Transition (A<->G or C<->T) vs transversion.
  [[nodiscard]] static bool isTransition(Alphabet nucleotide_1, Alphabet nucleotide_2) noexcept;

  /// Convert a base to an array offset.
  [[nodiscard]] static ContigOffset_t symbolToColumn(Alphabet nucleotide) noexcept;

  /// ReturnType nucleotide as a char.
  [[nodiscard]] static constexpr char convertToChar(Alphabet nucleotide) noexcept {
    return std::to_underlying(nucleotide);
  }

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CodingDNA5 — the stranded DNA policy. Shares the entire Nucleotide implementation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct CodingDNA5 {

  using Alphabet = Nucleotide;

  inline static constexpr ContigOffset_t NUCLEOTIDE_COLUMNS = 5;

  CodingDNA5() = delete;
  ~CodingDNA5() = delete;

  [[nodiscard]] static const std::array<Alphabet, 5>& enumerateAlphabet() noexcept;
  [[nodiscard]] static bool validAlphabet(Alphabet nucleotide) noexcept;
  [[nodiscard]] static ContigOffset_t symbolToColumn(Alphabet nucleotide) noexcept;
  [[nodiscard]] static Alphabet convertChar(char chr_base) { return DNA5::convertChar(chr_base); }

  [[nodiscard]] static constexpr char convertToChar(Alphabet nucleotide) noexcept {
    return std::to_underlying(nucleotide);
  }

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AminoAcid — the protein policy.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct AminoAcid {

  using Alphabet = Amino;

  inline static constexpr char START_CODON = 'M';
  inline static constexpr char STOP_CODON = '*';
  inline static constexpr char UNKNOWN_AMINO = 'Z';

  inline static constexpr Alphabet AMINO_STOP = Alphabet::Stop;
  inline static constexpr Alphabet AMINO_UNKNOWN = Alphabet::Z;

  AminoAcid() = delete;
  ~AminoAcid() = delete;

  /// ReturnType a vector of all valid alphabet values (20 natural + stop + unknown).
  [[nodiscard]] static const std::array<Alphabet, 22>& enumerateAlphabet() noexcept;

  /// Checks for possible memory corruption (includes the rare U/O amino acids).
  [[nodiscard]] static bool validAlphabet(Alphabet amino) noexcept;

  /// Convert an amino into an array offset.
  [[nodiscard]] static ContigOffset_t symbolToColumn(Alphabet amino) noexcept;

  /// Convert char to an Amino symbol (logs on invalid input).
  [[nodiscard]] static Alphabet convertChar(char chr_aa);

  /// ReturnType amino acid as a char.
  [[nodiscard]] static constexpr char convertToChar(Alphabet amino) noexcept {
    return std::to_underlying(amino);
  }

};


}   // end namespace


#endif //KGL_ALPHABET_H
