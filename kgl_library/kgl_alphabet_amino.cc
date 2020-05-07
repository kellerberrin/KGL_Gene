//
// Created by kellerberrin on 16/11/17.
//

#include "kgl_alphabet_amino.h"

namespace kgl = kellerberrin::genome;



// Covert from char to alphabet
kgl::AminoAcid::Alphabet kgl::AminoAcid::convertChar(char chr_base) {

  // Translate the nucleotide to an array column
  switch (chr_base) {

    case PHENYLALANINE: return Alphabet::F;
    case LEUCINE: return Alphabet::L;
    case SERINE: return Alphabet::S;
    case TYROSINE: return Alphabet::Y;
    case CYSTEINE: return Alphabet::C;
    case TRYPTOPHAN: return Alphabet::W;
    case PROLINE: return Alphabet::P;
    case HISTIDINE: return Alphabet::H;
    case GLUTAMINE: return Alphabet::Q;
    case ARGININE: return Alphabet::R;
    case ISOLEUCINE: return Alphabet::I;
    case METHIONINE: return Alphabet::M;
    case THREONINE: return Alphabet::T;
    case ASPARAGINE: return Alphabet::N;
    case LYSINE: return Alphabet::K;
    case VALINE: return Alphabet::V;
    case ALANINE: return Alphabet::A;
    case ASPARTIC: return Alphabet::D;
    case GLUTAMIC: return Alphabet::E;
    case GLYCINE: return Alphabet::G;
      // Rare - The additional two amino acids encoded using stop codons by some species.
    case SELENOCYSTEINE: return Alphabet::U;
    case PYRROLYSINE: return Alphabet::O;
      // The three stop codons.
    case STOP_CODON: return Alphabet::_;
      // The special unknown amino acid generated when the DNA5 codon
      // contains the unknown base 'N'.
    case UNKNOWN_AMINO: return Alphabet::Z;

    default:
      ExecEnv::log().error("AminoAcid::Alphabet(), Invalid nucleotide: {}", chr_base);
      return Alphabet::Z;

  }

}


// Covert alphabet symbol to an offset, used with the function above to create and access vectors of AA symbols.
kgl::ContigOffset_t kgl::AminoAcid::symbolToColumn(Alphabet amino) {

  // Translate the nucleotide to an array column
  switch (amino) {

    case Alphabet::F: return PHENYLALANINE_OFFSET;
    case Alphabet::L: return LEUCINE_OFFSET;
    case Alphabet::S: return SERINE_OFFSET;
    case Alphabet::Y: return TYROSINE_OFFSET;
    case Alphabet::C: return CYSTEINE_OFFSET;
    case Alphabet::W: return TRYPTOPHAN_OFFSET;
    case Alphabet::P: return PROLINE_OFFSET;
    case Alphabet::H: return HISTIDINE_OFFSET;
    case Alphabet::Q: return GLUTAMINE_OFFSET;
    case Alphabet::R: return ARGININE_OFFSET;
    case Alphabet::I: return ISOLEUCINE_OFFSET;
    case Alphabet::M: return METHIONINE_OFFSET;
    case Alphabet::T: return THREONINE_OFFSET;
    case Alphabet::N: return ASPARAGINE_OFFSET;
    case Alphabet::K: return LYSINE_OFFSET;
    case Alphabet::V: return VALINE_OFFSET;
    case Alphabet::A: return ALANINE_OFFSET;
    case Alphabet::D: return ASPARTIC_OFFSET;
    case Alphabet::E: return GLUTAMIC_OFFSET;
    case Alphabet::G: return GLYCINE_OFFSET;
      // The three stop codons.
    case Alphabet::_: return STOP_CODON_OFFSET;
      // The special unknown amino acid generated when the DNA5 codon
      // contains the unknown base 'N'.
    case Alphabet::Z: return UNKNOWN_AMINO_OFFSET;

    default:
      ExecEnv::log().error("AminoAcid::Alphabet(), Invalid amino symbol: {}", static_cast<char>(amino));
      return UNKNOWN_AMINO_OFFSET;

  }

}



bool kgl::AminoAcid::validAlphabet(Alphabet amino) {

  auto int_value = static_cast<size_t>(amino);

// We DO NOT use a switch here.
// Because the switch always assumes we can only have the enum values (and a memory corrupted sequence may not).
  bool compare = int_value == static_cast<size_t>(Alphabet::F)
                 or int_value == static_cast<size_t>(Alphabet::L)
                 or int_value == static_cast<size_t>(Alphabet::S)
                 or int_value == static_cast<size_t>(Alphabet::Y)
                 or int_value == static_cast<size_t>(Alphabet::C)
                 or int_value == static_cast<size_t>(Alphabet::W)
                 or int_value == static_cast<size_t>(Alphabet::P)
                 or int_value == static_cast<size_t>(Alphabet::H)
                 or int_value == static_cast<size_t>(Alphabet::Q)
                 or int_value == static_cast<size_t>(Alphabet::R)
                 or int_value == static_cast<size_t>(Alphabet::I)
                 or int_value == static_cast<size_t>(Alphabet::M)
                 or int_value == static_cast<size_t>(Alphabet::T)
                 or int_value == static_cast<size_t>(Alphabet::N)
                 or int_value == static_cast<size_t>(Alphabet::K)
                 or int_value == static_cast<size_t>(Alphabet::V)
                 or int_value == static_cast<size_t>(Alphabet::A)
                 or int_value == static_cast<size_t>(Alphabet::D)
                 or int_value == static_cast<size_t>(Alphabet::E)
                 or int_value == static_cast<size_t>(Alphabet::G)
//                 or int_value == static_cast<size_t>(Alphabet::U)
//                 or int_value == static_cast<size_t>(Alphabet::O)
                 or int_value == static_cast<size_t>(Alphabet::_)
                 or int_value == static_cast<size_t>(Alphabet::Z);

  return compare;

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
