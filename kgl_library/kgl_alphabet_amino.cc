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
                 or int_value == static_cast<size_t>(Alphabet::U)
                 or int_value == static_cast<size_t>(Alphabet::O)
                 or int_value == static_cast<size_t>(Alphabet::_)
                 or int_value == static_cast<size_t>(Alphabet::Z);

  return compare;

}


std::vector<kgl::AminoAcid::Alphabet> kgl::AminoAcid::enumerateAlphabet() {

  std::vector<Alphabet> alphabet_vector;

  alphabet_vector.emplace_back(Alphabet::F);
  alphabet_vector.emplace_back(Alphabet::L);
  alphabet_vector.emplace_back(Alphabet::S);
  alphabet_vector.emplace_back(Alphabet::Y);
  alphabet_vector.emplace_back(Alphabet::C);
  alphabet_vector.emplace_back(Alphabet::W);
  alphabet_vector.emplace_back(Alphabet::P);
  alphabet_vector.emplace_back(Alphabet::H);
  alphabet_vector.emplace_back(Alphabet::Q);
  alphabet_vector.emplace_back(Alphabet::R);
  alphabet_vector.emplace_back(Alphabet::I);
  alphabet_vector.emplace_back(Alphabet::M);
  alphabet_vector.emplace_back(Alphabet::T);
  alphabet_vector.emplace_back(Alphabet::N);
  alphabet_vector.emplace_back(Alphabet::K);
  alphabet_vector.emplace_back(Alphabet::V);
  alphabet_vector.emplace_back(Alphabet::A);
  alphabet_vector.emplace_back(Alphabet::D);
  alphabet_vector.emplace_back(Alphabet::E);
  alphabet_vector.emplace_back(Alphabet::G);
  // The special unknown amino acid generated when the DNA5 codon
  // contains the unknown base 'N'.
  alphabet_vector.emplace_back(Alphabet::Z);
  // Generic stop codon.
  alphabet_vector.emplace_back(Alphabet::_);

  return alphabet_vector;



}
