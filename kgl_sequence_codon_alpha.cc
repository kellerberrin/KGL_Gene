//
// Created by kellerberrin on 16/11/17.
//


#include "kgl_sequence_codon_alpha.h"


namespace kgl = kellerberrin::genome;



kgl::Codon::Codon(std::shared_ptr<DNA5SequenceCoding> sequence_ptr, ContigOffset_t codon_index) {

  if (codon_index >= codonLength(sequence_ptr)) {

    ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                         codon_index, sequence_ptr->length());
    codon_index = 0;
  }

  codon_index = static_cast<ContigOffset_t>(codon_index * 3);

  bases_[0] = sequence_ptr->at(codon_index);
  ++codon_index;
  bases_[1] = sequence_ptr->at(codon_index);
  ++codon_index;
  bases_[2] = sequence_ptr->at(codon_index);

}

