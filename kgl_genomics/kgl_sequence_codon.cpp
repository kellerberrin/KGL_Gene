//
// Created by kellerberrin on 14/11/17.
//


#include "kgl_sequence_codon.h"


namespace kgl = kellerberrin::genome;



kgl::Codon::Codon(const DNA5SequenceCoding& coding_sequence, ContigOffset_t codon_index) {

  if (codon_index >= codonLength(coding_sequence)) {

    ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                         codon_index, coding_sequence.length());
    codon_index = 0;

  }

  codon_index = static_cast<ContigOffset_t>(codon_index * Codon::CODON_SIZE);

  bases_[0] = coding_sequence.at(codon_index);
  ++codon_index;
  bases_[1] = coding_sequence.at(codon_index);
  ++codon_index;
  bases_[2] = coding_sequence.at(codon_index);

}


void kgl::Codon::modifyBase(size_t index, CodingDNA5::Alphabet base)
{

  if (index >= CODON_SIZE) {

    ExecEnv::log().error("Invalid codon base index specified index:{}, must be < 3)", index);
    return;
  }
  bases_[index] = base;

}
