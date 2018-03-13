//
// Created by kellerberrin on 15/02/18.
//



#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <edlib.h>


#include "kgl_sequence_compare_impl.h"
#include "kgl_exec_env.h"
#include "kgl_genome_types.h"
#include "kgl_genome_db.h"
#include "kgl_sequence_base.h"
#include "kgl_sequence_codon.h"


namespace kgl = kellerberrin::genome;


class kgl::SequenceComparison::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;


  kgl::CompareScore_t MyerHirschbergGlobal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;
  kgl::CompareScore_t MyerHirschbergLocal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;
  kgl::CompareScore_t DNALocalAffineGap(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  std::string editItems(const std::string& reference_str, const std::string& mutant_str, char delimiter, VariantOutputIndex display_offset) const;

  void editDNAItems(std::shared_ptr<const ContigFeatures> contig_ptr,
                    std::shared_ptr<const DNA5SequenceCoding> reference,
                    std::shared_ptr<const DNA5SequenceCoding> mutant,
                    MutationEditVector& mutation_edit_vector) const;
private:

  EditVector createEditItems(const std::string& reference_str, const std::string& mutant_str, EditVector& edit_vector) const;
  EditVector createEditItemsEdlib(const std::string& reference_str, const std::string& mutant_str, EditVector& edit_vector) const;

};





kgl::CompareScore_t kgl::SequenceComparison::SequenceManipImpl::MyerHirschbergGlobal(const std::string& sequenceA,
                                                                                       const std::string& sequenceB,
                                                                                       std::string& compare_str) const {

  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;

  TSequence seq1 = sequenceA.c_str();
  TSequence seq2 = sequenceB.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  CompareScore_t score = seqan::globalAlignment(align, seqan::MyersHirschberg());
  ss << align << std::endl;

  compare_str = ss.str();

  return score;

}



kgl::CompareScore_t kgl::SequenceComparison::SequenceManipImpl::MyerHirschbergLocal(const std::string& sequenceA,
                                                                                      const std::string& sequenceB,
                                                                                      std::string& compare_str) const {

  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;

  TSequence seq1 = sequenceA.c_str();
  TSequence seq2 = sequenceB.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  CompareScore_t score = seqan::globalAlignment(align, seqan::MyersHirschberg());
  ss << align << std::endl;

  compare_str = ss.str();

  return score;

}




kgl::CompareScore_t kgl::SequenceComparison::SequenceManipImpl::DNALocalAffineGap(const std::string& sequenceA,
                                                                                    const std::string& sequenceB,
                                                                                    std::string& compare_str) const {

  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;

  TSequence seq1 = sequenceA.c_str();
  TSequence seq2 = sequenceB.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

//  CompareScore_t score = seqan::localAlignment(align, seqan::Score<int>(0, -1, -2, -2), seqan::DynamicGaps());
//  CompareScore_t score = seqan::globalAlignment(align, seqan::Score<int>(0, -1, -2, -1), seqan::AlignConfig<true, false, false, true>(), seqan::AffineGaps());
//  CompareScore_t score = seqan::globalAlignment(align, seqan::Score<int>(3, -3, -20, -20), seqan::AlignConfig<true, false, false, true>(), seqan::AffineGaps());
  CompareScore_t score = seqan::globalAlignment(align, seqan::Score<int>(3, -3, -5, -2), seqan::AlignConfig<true, false, false, true>(), seqan::AffineGaps());


  ss << align << std::endl;

  compare_str = ss.str();

  return score;

}


void kgl::SequenceComparison::SequenceManipImpl::editDNAItems(std::shared_ptr<const ContigFeatures> contig_ptr,
                                                              std::shared_ptr<const DNA5SequenceCoding> reference,
                                                              std::shared_ptr<const DNA5SequenceCoding> mutant,
                                                              MutationEditVector& mutation_edit_vector) const {
  EditVector edit_vector;
  EditVector check_edit_vector;
  createEditItems(reference->getSequenceAsString(), mutant->getSequenceAsString(), edit_vector);
  createEditItemsEdlib(reference->getSequenceAsString(), mutant->getSequenceAsString(), check_edit_vector);

  bool match = true;
  if (edit_vector.size() != check_edit_vector.size()) {

    ExecEnv::log().error("Edit vector size mis-match, edit_vector size: {}, check edit vector size: {}", edit_vector.size(), check_edit_vector.size());
    match = false;

  } else {

    for (size_t index = 0; index < edit_vector.size(); ++index) {

      bool compare = edit_vector[index].mutant_char == check_edit_vector[index].mutant_char
                     and edit_vector[index].reference_offset == check_edit_vector[index].reference_offset
                     and edit_vector[index].reference_char == check_edit_vector[index].reference_char;


      if (not compare) {

        ExecEnv::log().error("Edit item mis-match, edit_vector: {}{}{}, check edit vector: {}{}{}",
                             edit_vector[index].reference_char, edit_vector[index].reference_offset, edit_vector[index].mutant_char,
                             check_edit_vector[index].reference_char, check_edit_vector[index].reference_offset, check_edit_vector[index].mutant_char);
        match = false;
      }

    }

  }

  if (not match) {

    std::string compare_str;
    MyerHirschbergGlobal(reference->getSequenceAsString(), mutant->getSequenceAsString(), compare_str);
    ExecEnv::log().info("Edit vector mis-match sequence comparison:\n {}", compare_str);

  }

  mutation_edit_vector.mutation_vector.clear();
  MutationItem mutation_item;
  for (auto edit_item :edit_vector) {

    ContigOffset_t codon_index = static_cast<size_t>(edit_item.reference_offset / 3);

    std::shared_ptr<const Codon> mutant_codon(std::make_shared<Codon>(mutant, codon_index));
    std::shared_ptr<const Codon> ref_codon(std::make_shared<Codon>(reference, codon_index));

    mutation_item.DNA_mutation = edit_item;
    mutation_item.reference_codon = ref_codon->getSequenceAsString();
    mutation_item.mutation_codon = mutant_codon->getSequenceAsString();

    mutation_item.amino_mutation.reference_char = AminoAcid::convertToChar(contig_ptr->getAminoAcid(*ref_codon));
    mutation_item.amino_mutation.reference_offset = codon_index;
    mutation_item.amino_mutation.mutant_char = AminoAcid::convertToChar(contig_ptr->getAminoAcid(*mutant_codon));

    mutation_edit_vector.mutation_vector.push_back(mutation_item);

  }

}



std::string kgl::SequenceComparison::SequenceManipImpl::editItems(const std::string& reference,
                                                                    const std::string& mutant,
                                                                    char delimiter,
                                                                    VariantOutputIndex index_offset) const {
  EditVector edit_vector;
  createEditItems(reference, mutant, edit_vector);
  size_t last_index = edit_vector.size() - 1;
  size_t index = 0;
  std::stringstream ss;
  for (auto edit_item :edit_vector) {

    ss << edit_item.reference_char << offsetOutput(edit_item.reference_offset, index_offset) << edit_item.mutant_char;

    if (index != last_index) {

      ss << delimiter;

    }

    index++;


  }

  return ss.str();

}


kgl::EditVector kgl::SequenceComparison::SequenceManipImpl::createEditItemsEdlib(const std::string& reference,
                                                                            const std::string& mutant,
                                                                            EditVector& edit_vector) const {

  edit_vector.clear();

  EdlibAlignResult result = edlibAlign(mutant.c_str(), mutant.size(),reference.c_str(), reference.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
  if (result.status == EDLIB_STATUS_OK) {

    EditItem edit_item;
    int mut_index = 0;
    for (int ref_index = 0; ref_index < result.alignmentLength; ++ref_index) {

      switch(result.alignment[ref_index]) {

        case 0:
          break;

        case 1:
          edit_item.reference_char = reference[ref_index];
          edit_item.reference_offset = ref_index;
          edit_item.mutant_char = EditItem::INSERTION;
          edit_vector.push_back(edit_item);
          ++mut_index;
          break;

        case 2:
          edit_item.reference_char = reference[ref_index];
          edit_item.reference_offset = ref_index;
          edit_item.mutant_char = EditItem::DELETION;
          edit_vector.push_back(edit_item);
          --mut_index;
          break;

        case 3:
          edit_item.reference_char = reference[ref_index];
          edit_item.reference_offset = ref_index;
          edit_item.mutant_char = mutant[mut_index];
          edit_vector.push_back(edit_item);
          break;

      }

      ++mut_index;

    }


  } else {

    ExecEnv::log().error("Edlib - problem generating cigar reference:{}, alternate: {}", reference, mutant);

  }

  edlibFreeAlignResult(result);

  return edit_vector;

}



kgl::EditVector kgl::SequenceComparison::SequenceManipImpl::createEditItems(const std::string& reference_str,
                                                                         const std::string& mutant_str,
                                                                         EditVector& edit_vector) const
{
  typedef seqan::String<char> TSequence;
  typedef seqan::Gaps<TSequence, seqan::ArrayGaps> TGaps;
  typedef seqan::Iterator<TGaps>::Type TGapsIterator;


  TSequence reference = reference_str.c_str();
  TSequence mutant = mutant_str.c_str();

  TGaps gaps_reference;
  TGaps gaps_mutant;

  seqan::assignSource(gaps_reference, reference);
  seqan::assignSource(gaps_mutant, mutant);

  seqan::clearGaps(gaps_reference);
  seqan::clearGaps(gaps_mutant);

  seqan::globalAlignment(gaps_reference, gaps_mutant, seqan::Score<int>(1, -1, -10, -7));

  TGapsIterator itGapsText = begin(gaps_mutant);
  TGapsIterator itGapsPattern = begin(gaps_reference);
  TGapsIterator itGapsEnd = end(gaps_reference);


  ContigOffset_t ref_index = 0;
  ContigOffset_t mut_index = 0;

  while (itGapsPattern != itGapsEnd)
  {

    // Count insertions.
    if (seqan::isGap(itGapsText))
    {
      int numGaps = seqan::countGaps(itGapsText);
      itGapsText += numGaps;
      itGapsPattern += numGaps;

      for(int i = 0; i < numGaps; ++i) {

        EditItem edit_item;
        edit_item.reference_char = reference_str[ref_index];
        edit_item.reference_offset = ref_index;
        edit_item.mutant_char = EditItem::DELETION;
        edit_vector.push_back(edit_item);
        ++ref_index;

      }
      continue;
    }

    // Count deletions.
    if (seqan::isGap(itGapsPattern))
    {
      int numGaps = seqan::countGaps(itGapsPattern);
      itGapsText += numGaps;
      itGapsPattern += numGaps;
      for(int i = 0; i < numGaps; ++i) {

        EditItem edit_item;
        edit_item.reference_char = reference_str[ref_index];
        edit_item.reference_offset = ref_index;
        edit_item.mutant_char = EditItem::INSERTION;
        edit_vector.push_back(edit_item);
        ++mut_index;

      }
      continue;
    }

    // Count matches and  mismatches.
    while (itGapsPattern != itGapsEnd)
    {
      if (seqan::isGap(itGapsPattern) || seqan::isGap(itGapsText))
        break;

      if (reference_str[ref_index] != mutant_str[mut_index]) {

        EditItem edit_item;
        edit_item.reference_char = reference_str[ref_index];
        edit_item.reference_offset = ref_index;
        edit_item.mutant_char = mutant_str[mut_index];
        edit_vector.push_back(edit_item);

      }
      ++ref_index;
      ++mut_index;
      ++itGapsText;
      ++itGapsPattern;
    }

  }
  // Output the hit position in the text, the total number of edits and the corresponding cigar string.

  return edit_vector;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::FreeBayesVCFImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::SequenceComparison::SequenceComparison() : sequence_manip_impl_ptr_(std::make_unique<kgl::SequenceComparison::SequenceManipImpl>()) {}
kgl::SequenceComparison::~SequenceComparison() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.

// Comparison

kgl::CompareScore_t kgl::SequenceComparison::MyerHirschbergGlobal(const std::string& sequenceA,
                                                                    const std::string& sequenceB,
                                                                    std::string& compare_str) const {

  return sequence_manip_impl_ptr_->MyerHirschbergGlobal(sequenceA, sequenceB, compare_str);

}


kgl::CompareScore_t kgl::SequenceComparison::MyerHirschbergLocal(const std::string& sequenceA,
                                                                    const std::string& sequenceB,
                                                                    std::string& compare_str) const {

  return sequence_manip_impl_ptr_->MyerHirschbergLocal(sequenceA, sequenceB, compare_str);

}



kgl::CompareScore_t kgl::SequenceComparison::DNALocalAffineGap(const std::string& sequenceA,
                                                                   const std::string& sequenceB,
                                                                   std::string& compare_str) const {

  return sequence_manip_impl_ptr_->DNALocalAffineGap(sequenceA, sequenceB, compare_str);

}


std::string kgl::SequenceComparison::editItems(const std::string& reference,
                                                 const std::string& mutant,
                                                 char delimiter,
                                                 VariantOutputIndex index_offset) const {

  return sequence_manip_impl_ptr_->editItems(reference, mutant, delimiter, index_offset);

}



void kgl::SequenceComparison::editDNAItems(std::shared_ptr<const ContigFeatures> contig_ptr,
                                           std::shared_ptr<const DNA5SequenceCoding> reference,
                                           std::shared_ptr<const DNA5SequenceCoding> mutant,
                                           MutationEditVector& mutation_edit_vector) const {

  sequence_manip_impl_ptr_->editDNAItems(contig_ptr, reference, mutant, mutation_edit_vector);

}



std::string kgl::SequenceComparison::editDNAItems(std::shared_ptr<const ContigFeatures> contig_ptr,
                                                  std::shared_ptr<const DNA5SequenceCoding> reference,
                                                  std::shared_ptr<const DNA5SequenceCoding> mutant,
                                                  char delimiter,
                                                  VariantOutputIndex index_offset) const {


  MutationEditVector mutation_edit_vector;
  editDNAItems(contig_ptr, reference, mutant, mutation_edit_vector);

  size_t last_index = mutation_edit_vector.mutation_vector.size() - 1;
  size_t index = 0;
  std::stringstream ss;
  for (auto edit_item : mutation_edit_vector.mutation_vector) {

    ss << edit_item.DNA_mutation.reference_char
       << offsetOutput(edit_item.DNA_mutation.reference_offset, index_offset)
       << edit_item.DNA_mutation.mutant_char;
    ss << " " << edit_item.reference_codon << "-" << edit_item.mutation_codon;
    ss << " " << edit_item.amino_mutation.reference_char
       << offsetOutput(edit_item.amino_mutation.reference_offset, index_offset)
       << edit_item.amino_mutation.mutant_char;

    if (index != last_index) {

      ss << delimiter;

    }

    index++;


  }

  return ss.str();

}


