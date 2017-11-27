//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A compound variant. A collection of feature aligned and contiguous variants. Insertions and Deletions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::CompoundVariant::equivalent(const Variant& cmp_var) const {

  auto compound_var = dynamic_cast<const CompoundVariant*>(&cmp_var);

  if (compound_var == nullptr) return false;

  bool result = contigId() == compound_var->contigId()
                and contigOffset() == compound_var->contigOffset()
                and type() == compound_var->type()
                and codingSequenceId() == compound_var->codingSequenceId();

  if (not result) return false;

  if (getMap().size() != compound_var->getMap().size()) return false;

  auto cmp_iterator = compound_var->getMap().begin();
  for (auto iterator = getMap().begin(); iterator != getMap().end(); ++iterator) {

    if(cmp_iterator == compound_var->getMap().end()) return false;

    if (not iterator->second->equivalent(*(cmp_iterator->second))) return false;

    cmp_iterator++;

  }

  return true;

}


std::string kgl::CompoundVariant::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << variant_map_.size() << delimiter;
  ss << location(delimiter, output_index);
  ss << mutation(delimiter, output_index) << "\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->suboutput(delimiter, output_index);
  }
  return ss.str();

}

std::string kgl::CompoundVariant::location(char delimiter, VariantOutputIndex output_index) const {

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

      ExecEnv::log().error("Compound variant in contig: {} offset: {} is empty", contigId(), offset());

    }

  } else {

    ExecEnv::log().error("Compound variant in contig: {} offset: {} is not in a coding sequence", contigId(), offset());

  }

  return ss.str();

}

