//
// Created by kellerberrin on 24/11/17.
//


#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;

std::string kgl::CompoundInsert::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << variant_map_.size() << delimiter;
  ss << mutation(delimiter, output_index) << "\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->suboutput(delimiter, output_index);
  }
  return ss.str();

}

std::string kgl::CompoundInsert::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    if (not getMap().empty()) {

      ContigSize_t base_in_codon_begin;
      ContigOffset_t codon_offset_begin;
      getMap().begin()->second->codonOffset(codon_offset_begin, base_in_codon_begin);
      ContigSize_t base_in_codon_end;
      ContigOffset_t codon_offset_end;
      getMap().rbegin()->second->codonOffset(codon_offset_end, base_in_codon_end);

      if (codon_offset_begin == codon_offset_end) {

        if (base_in_codon_begin < base_in_codon_end) {

          ss << offsetOutput(codon_offset_begin, output_index) << CODON_BASE_SEPARATOR
             << offsetOutput(base_in_codon_begin, output_index) << delimiter;

        } else {

          ss << offsetOutput(codon_offset_end, output_index) << CODON_BASE_SEPARATOR
             << offsetOutput(base_in_codon_end, output_index) << delimiter;

        }

      } else if (codon_offset_begin < codon_offset_end) {

        ss << offsetOutput(codon_offset_begin, output_index) << CODON_BASE_SEPARATOR
           << offsetOutput(base_in_codon_begin, output_index) << delimiter;

      } else {

        ss << offsetOutput(codon_offset_end, output_index) << CODON_BASE_SEPARATOR
           << offsetOutput(base_in_codon_end, output_index) << delimiter;

      }


    } else {

      ExecEnv::log().error("Compound insert in contig: {} offset: {} is empty", contigId(), offset());

    }

  } else {

    ExecEnv::log().error("Compound insert in contig: {} offset: {} is not in a coding sequence", contigId(), offset());

  }

  return ss.str();

}


bool kgl::CompoundInsert::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                               std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {

  return true;

}

