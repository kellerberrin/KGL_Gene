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

  ContigOffset_t base_offset = codon_index * CODON_SIZE;

  bases_[0] = coding_sequence.at(base_offset);
  bases_[1] = coding_sequence.at(base_offset + 1);
  bases_[2] = coding_sequence.at(base_offset + 2);

}