//
// Created by kellerberrin on 16/11/17.
//

#include "kgl_alphabet_amino.h"
#include "kel_exec_env.h"

namespace kgl = kellerberrin::genome;


namespace kellerberrin::genome::detail {   //  Shared constexpr tables for the amino alphabet.


// Total lookup tables over the unsigned char domain: branch-free and a memory-corrupted
// raw value can never index out of bounds (the no-switch intent of the original code).


inline constexpr auto AMINO_COLUMN_TABLE = []() consteval {

  std::array<ContigOffset_t, 256> table{};
  table.fill(AminoAcid::UNKNOWN_AMINO_OFFSET);
  table[static_cast<unsigned char>(AminoAcid::PHENYLALANINE)] = AminoAcid::PHENYLALANINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::LEUCINE)] = AminoAcid::LEUCINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::SERINE)] = AminoAcid::SERINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::TYROSINE)] = AminoAcid::TYROSINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::CYSTEINE)] = AminoAcid::CYSTEINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::TRYPTOPHAN)] = AminoAcid::TRYPTOPHAN_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::PROLINE)] = AminoAcid::PROLINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::HISTIDINE)] = AminoAcid::HISTIDINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::GLUTAMINE)] = AminoAcid::GLUTAMINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::ARGININE)] = AminoAcid::ARGININE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::ISOLEUCINE)] = AminoAcid::ISOLEUCINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::METHIONINE)] = AminoAcid::METHIONINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::THREONINE)] = AminoAcid::THREONINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::ASPARAGINE)] = AminoAcid::ASPARAGINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::LYSINE)] = AminoAcid::LYSINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::VALINE)] = AminoAcid::VALINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::ALANINE)] = AminoAcid::ALANINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::ASPARTIC)] = AminoAcid::ASPARTIC_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::GLUTAMIC)] = AminoAcid::GLUTAMIC_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::GLYCINE)] = AminoAcid::GLYCINE_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::STOP_CODON)] = AminoAcid::STOP_CODON_OFFSET;
  table[static_cast<unsigned char>(AminoAcid::UNKNOWN_AMINO)] = AminoAcid::UNKNOWN_AMINO_OFFSET;
  return table;

}();

// B3a: includes the rare selenocysteine (U) and pyrrolysine (O) so a parsed sequence
// containing them passes the module's own verifyString() corruption check.
inline constexpr auto AMINO_VALID_TABLE = []() consteval {

  std::array<bool, 256> table{};
  for (size_t index = 0; index < 256; ++index) {

    table[index] = (index == static_cast<unsigned char>(AminoAcid::UNKNOWN_AMINO))
                   or (AMINO_COLUMN_TABLE[index] != AminoAcid::UNKNOWN_AMINO_OFFSET);

  }
  table[static_cast<unsigned char>(AminoAcid::SELENOCYSTEINE)] = true;
  table[static_cast<unsigned char>(AminoAcid::PYRROLYSINE)] = true;
  return table;

}();

inline constexpr auto AMINO_ALPHABET_TABLE = []() consteval {

  using Alphabet = AminoAcid::Alphabet;
  std::array<Alphabet, 256> table{};
  table.fill(Alphabet::Z);
  table[static_cast<unsigned char>(AminoAcid::PHENYLALANINE)] = Alphabet::F;
  table[static_cast<unsigned char>(AminoAcid::LEUCINE)] = Alphabet::L;
  table[static_cast<unsigned char>(AminoAcid::SERINE)] = Alphabet::S;
  table[static_cast<unsigned char>(AminoAcid::TYROSINE)] = Alphabet::Y;
  table[static_cast<unsigned char>(AminoAcid::CYSTEINE)] = Alphabet::C;
  table[static_cast<unsigned char>(AminoAcid::TRYPTOPHAN)] = Alphabet::W;
  table[static_cast<unsigned char>(AminoAcid::PROLINE)] = Alphabet::P;
  table[static_cast<unsigned char>(AminoAcid::HISTIDINE)] = Alphabet::H;
  table[static_cast<unsigned char>(AminoAcid::GLUTAMINE)] = Alphabet::Q;
  table[static_cast<unsigned char>(AminoAcid::ARGININE)] = Alphabet::R;
  table[static_cast<unsigned char>(AminoAcid::ISOLEUCINE)] = Alphabet::I;
  table[static_cast<unsigned char>(AminoAcid::METHIONINE)] = Alphabet::M;
  table[static_cast<unsigned char>(AminoAcid::THREONINE)] = Alphabet::T;
  table[static_cast<unsigned char>(AminoAcid::ASPARAGINE)] = Alphabet::N;
  table[static_cast<unsigned char>(AminoAcid::LYSINE)] = Alphabet::K;
  table[static_cast<unsigned char>(AminoAcid::VALINE)] = Alphabet::V;
  table[static_cast<unsigned char>(AminoAcid::ALANINE)] = Alphabet::A;
  table[static_cast<unsigned char>(AminoAcid::ASPARTIC)] = Alphabet::D;
  table[static_cast<unsigned char>(AminoAcid::GLUTAMIC)] = Alphabet::E;
  table[static_cast<unsigned char>(AminoAcid::GLYCINE)] = Alphabet::G;
  table[static_cast<unsigned char>(AminoAcid::SELENOCYSTEINE)] = Alphabet::U;
  table[static_cast<unsigned char>(AminoAcid::PYRROLYSINE)] = Alphabet::O;
  table[static_cast<unsigned char>(AminoAcid::STOP_CODON)] = Alphabet::_;
  table[static_cast<unsigned char>(AminoAcid::UNKNOWN_AMINO)] = Alphabet::Z;
  return table;

}();


}   // end namespace


// Covert from char to alphabet.
kgl::AminoAcid::Alphabet kgl::AminoAcid::convertChar(char chr_base) {

  const Alphabet amino = detail::AMINO_ALPHABET_TABLE[static_cast<unsigned char>(chr_base)];
  if (amino == Alphabet::Z and chr_base != UNKNOWN_AMINO) {

    ExecEnv::log().error("AminoAcid::convertChar(), Invalid amino acid: '{}', ascii value: {}. Input is probably corrupt or not protein text.", chr_base, static_cast<size_t>(chr_base));

  }

  return amino;

}


// Covert alphabet symbol to an offset, used with the function above to create and access vectors of AA symbols.
kgl::ContigOffset_t kgl::AminoAcid::symbolToColumn(Alphabet amino) {

  // The table lookup is total so a corrupted raw value cannot index out of bounds.
  const ContigOffset_t column = detail::AMINO_COLUMN_TABLE[static_cast<unsigned char>(amino)];

  // The reference logged a corruption diagnostic here; restore it in table form.
  // validAlphabet() includes the B3a-approved U/O (silent); only corrupted values log.
  if (not validAlphabet(amino)) {

    ExecEnv::log().error("AminoAcid::symbolToColumn(), Invalid amino symbol: {}", static_cast<char>(amino));

  }

  return column;

}


bool kgl::AminoAcid::validAlphabet(Alphabet amino) {

  // We DO NOT use a switch here.
  // Because the switch always assumes we can only have the enum values (and a memory corrupted sequence may not).
  return detail::AMINO_VALID_TABLE[static_cast<unsigned char>(amino)];

}


const std::vector<kgl::AminoAcid::Alphabet>& kgl::AminoAcid::enumerateAlphabet() {

  static std::vector<Alphabet> alphabet_vector = { Alphabet::F,
                                                   Alphabet::L,
                                                   Alphabet::S,
                                                   Alphabet::Y,
                                                   Alphabet::C,
                                                   Alphabet::W,
                                                   Alphabet::P,
                                                   Alphabet::H,
                                                   Alphabet::Q,
                                                   Alphabet::R,
                                                   Alphabet::I,
                                                   Alphabet::M,
                                                   Alphabet::T,
                                                   Alphabet::N,
                                                   Alphabet::K,
                                                   Alphabet::V,
                                                   Alphabet::A,
                                                   Alphabet::D,
                                                   Alphabet::E,
                                                   Alphabet::G,
                                                   Alphabet::_, // Generic stop codon.
                                                   Alphabet::Z };   // Unknown amino acid generated when the DNA5 codon contains the unknown base 'N'.

  return alphabet_vector;

}