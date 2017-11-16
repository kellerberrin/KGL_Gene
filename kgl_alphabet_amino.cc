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
