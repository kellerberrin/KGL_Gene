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


namespace kgl = kellerberrin::genome;


struct EditItem{

  char reference_char;
  kgl::ContigOffset_t reference_offset;
  char mutant_char;

};

using EditVector = std::vector<EditItem>;


class kgl::SequenceManipulation::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;

  std::string blosum80(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  std::string localBlosum80(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  std::string localBlosum62(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score, ContigSize_t& length) const;

  std::string simplealign(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  kgl::CompareScore_t MyerHirschberg(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  std::string multipleAlign(const std::vector<std::string>& compare_str_vec) const;

  kgl::CompareScore_t MyerHirschbergGlobal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;
  kgl::CompareScore_t MyerHirschbergLocal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  kgl::CompareDistance_t LevenshteinGlobal(const std::string& sequenceA, const std::string& sequenceB) const;
  kgl::CompareDistance_t LevenshteinLocal(const std::string& sequenceA, const std::string& sequenceB) const;

  kgl::CompareDistance_t globalblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;


  std::string editItems(const std::string& reference_str, const std::string& mutant_str, char delimiter, VariantOutputIndex display_offset);


private:

  EditVector createEditItems(const std::string& reference_str, const std::string& mutant_str, EditVector& edit_vector);

};




std::string kgl::SequenceManipulation::SequenceManipImpl::localBlosum80(const std::string& reference_str,
                                                                        const std::string& compare_str,
                                                                        CompareScore_t& score) const {
  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;
  int gap = -1;

  TSequence seq1 = reference_str.c_str();
  TSequence seq2 = compare_str.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  score = seqan::localAlignment(align, seqan::Blosum80(gap, gap));
  ss << align << std::endl;

  return ss.str();

}




std::string kgl::SequenceManipulation::SequenceManipImpl::blosum80(const std::string& reference_str,
                                                                   const std::string& compare_str,
                                                                   CompareScore_t& score) const {
  using TSequence = seqan::String<seqan::AminoAcid> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;
  int open_gap = -8;
  int extend_gap = -3;

  TSequence seq1 = reference_str.c_str();
  TSequence seq2 = compare_str.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  score = seqan::globalAlignment(align, seqan::Blosum80(open_gap, extend_gap));
  ss << align << std::endl;

  return ss.str();

}




std::string kgl::SequenceManipulation::SequenceManipImpl::localBlosum62(const std::string& reference_str,
                                                                        const std::string& compare_str,
                                                                        CompareScore_t& score,
                                                                        ContigSize_t& length) const {
  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;
  int gap = -1;

  length = 1;

  TSequence seq1 = reference_str.c_str();
  TSequence seq2 = compare_str.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  score = seqan::localAlignment(align, seqan::Blosum62(gap, gap));
  ss << align << std::endl;

  return ss.str();

}


std::string kgl::SequenceManipulation::SequenceManipImpl::simplealign(const std::string& reference_str,
                                                                      const std::string& compare_str,
                                                                      CompareScore_t& score) const {
  int match = 0;
  int mismatch = -1;
  int gap = -1;

  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;

  TSequence seq1 = reference_str.c_str();
  TSequence seq2 = compare_str.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  score = seqan::globalAlignment(align, seqan::Score<int, seqan::Simple>(match, mismatch, gap));

  ss << align << std::endl;

  return ss.str();

}


std::string kgl::SequenceManipulation::SequenceManipImpl::multipleAlign(const std::vector<std::string>& compare_str_vec) const {

  using TSequence = seqan::String<char>;
//  using DNATAlign = seqan::Align <seqan::DnaString>;
//  using AminoAlign = seqan::Align<seqan::AminoAcid>;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;

  TAlign align;

  resize(rows(align), compare_str_vec.size());
  std::vector<TSequence> seq_array;

  for (size_t i = 0; i < compare_str_vec.size(); ++i) {

    TSequence seq1 = compare_str_vec[i].c_str();
    seq_array.push_back(seq1);
    seqan::assignSource(row(align, i), seq_array.back());

  }

  std::stringstream ss;
  seqan::globalMsaAlignment(align, seqan::SimpleScore(5, -3, -1, -3));
  ss << align << std::endl;

  return ss.str();

}




kgl::CompareDistance_t kgl::SequenceManipulation::SequenceManipImpl::globalblosum80Distance(const std::string& sequenceA,
                                                                                            const std::string& sequenceB) const {
  using TSequence = seqan::String<seqan::AminoAcid> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;
  int open_gap = -8;
  int extend_gap = -3;

  TSequence seq1 = sequenceA.c_str();
  TSequence seq2 = sequenceB.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  long score = seqan::globalAlignment(align, seqan::Blosum80(open_gap, extend_gap));

  return static_cast<double>(score) * -1.0;  // Invert the scores.

}


kgl::CompareScore_t kgl::SequenceManipulation::SequenceManipImpl::MyerHirschbergGlobal(const std::string& sequenceA,
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



kgl::CompareScore_t kgl::SequenceManipulation::SequenceManipImpl::MyerHirschbergLocal(const std::string& sequenceA,
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


kgl::CompareDistance_t kgl::SequenceManipulation::SequenceManipImpl::LevenshteinGlobal(const std::string& sequenceA,
                                                                                       const std::string& sequenceB) const {

  EdlibAlignResult result = edlibAlign(sequenceA.c_str(),
                                       sequenceA.size(),
                                       sequenceB.c_str(),
                                       sequenceB.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_LOC, NULL, 0));

  if (result.status != EDLIB_STATUS_OK) {

    kgl::ExecEnv::log().error("Problem calculating Global Levenshtein distance using edlib; sequenceA: {}, sequenceB: {}",
                              sequenceA, sequenceB);
    edlibFreeAlignResult(result);
    return 0;
  }

  kgl::CompareDistance_t distance = std::fabs(result.editDistance);
  edlibFreeAlignResult(result);
  return distance;

}


kgl::CompareDistance_t kgl::SequenceManipulation::SequenceManipImpl::LevenshteinLocal(const std::string& sequenceA,
                                                                                      const std::string& sequenceB) const {

  EdlibAlignResult result = edlibAlign(sequenceA.c_str(),
                                       sequenceA.size(),
                                       sequenceB.c_str(),
                                       sequenceB.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));

  if (result.status != EDLIB_STATUS_OK) {

    kgl::ExecEnv::log().error("Problem calculating Local Levenshtein distance using edlib; sequenceA: {}, sequenceB: {}",
                              sequenceA, sequenceB);
    edlibFreeAlignResult(result);
    return 0;
  }

  kgl::CompareDistance_t distance = std::fabs(result.editDistance);
  if (result.alignmentLength <= 0) {

    kgl::ExecEnv::log().error("Problem calculating Local Levenshtein alignment, length zero or -ve {}; sequenceA: {}, sequenceB: {}",
                              result.alignmentLength, sequenceA, sequenceB);
    edlibFreeAlignResult(result);
    return 0;

  }
  edlibFreeAlignResult(result);
  return distance / static_cast<double>(result.alignmentLength);


}


std::string kgl::SequenceManipulation::SequenceManipImpl::editItems(const std::string& reference,
                                                                    const std::string& mutant,
                                                                    char delimiter,
                                                                    VariantOutputIndex index_offset) {
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



EditVector kgl::SequenceManipulation::SequenceManipImpl::createEditItems(const std::string& reference_str,
                                                                         const std::string& mutant_str,
                                                                         EditVector& edit_vector)
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

  seqan::globalAlignment(gaps_reference, gaps_mutant, seqan::Score<int>(0, -1, -1));

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
        edit_item.mutant_char = '-';
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
        edit_item.mutant_char = '+';
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


kgl::SequenceManipulation::SequenceManipulation() : sequence_manip_impl_ptr_(std::make_unique<kgl::SequenceManipulation::SequenceManipImpl>()) {}
kgl::SequenceManipulation::~SequenceManipulation() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


kgl::CompareScore_t kgl::SequenceManipulation::MyerHirschbergGlobal(const std::string& sequenceA,
                                                                    const std::string& sequenceB,
                                                                    std::string& compare_str) const {

  return sequence_manip_impl_ptr_->MyerHirschbergGlobal(sequenceA, sequenceB, compare_str);

}


kgl::CompareScore_t kgl::SequenceManipulation::MyerHirschbergLocal(const std::string& sequenceA,
                                                                    const std::string& sequenceB,
                                                                    std::string& compare_str) const {

  return sequence_manip_impl_ptr_->MyerHirschbergLocal(sequenceA, sequenceB, compare_str);

}


kgl::CompareDistance_t kgl::SequenceManipulation::LevenshteinGlobal(const std::string& sequenceA,
                                                                    const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->LevenshteinGlobal(sequenceA, sequenceB);

}


kgl::CompareDistance_t kgl::SequenceManipulation::LevenshteinLocal(const std::string& sequenceA,
                                                                   const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->LevenshteinLocal(sequenceA, sequenceB);

}


kgl::CompareDistance_t kgl::SequenceManipulation::globalblosum80Distance(const std::string& sequenceA,
                                                                         const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->globalblosum80Distance(sequenceA, sequenceB);

}


std::string kgl::SequenceManipulation::editItems(const std::string& reference,
                                                 const std::string& mutant,
                                                 char delimiter,
                                                 VariantOutputIndex index_offset) {

  return sequence_manip_impl_ptr_->editItems(reference, mutant, delimiter, index_offset);

}
