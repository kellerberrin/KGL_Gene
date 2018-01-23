//
// Created by kellerberrin on 4/01/18.
//

#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include "kgl_sequence_manip.h"


namespace kgl = kellerberrin::genome;



class kgl::SequenceManipulation::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;

  std::string blosum80(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  std::string simplealign(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  std::string myersHirschberg(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  std::string multipleAlign(const std::vector<std::string>& compare_str_vec) const;

private:


};




std::string kgl::SequenceManipulation::SequenceManipImpl::blosum80(const std::string& reference_str,
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

  score = seqan::globalAlignment(align, seqan::Blosum80(gap, gap));
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



std::string kgl::SequenceManipulation::SequenceManipImpl::myersHirschberg(const std::string& reference_str,
                                                                          const std::string& compare_str,
                                                                          CompareScore_t& score) const {

  using TSequence = seqan::String<char> ;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps> ;

  TSequence seq1 = reference_str.c_str();
  TSequence seq2 = compare_str.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan::assignSource(row(align, 0), seq1);
  seqan::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  score = seqan::globalAlignment(align, seqan::MyersHirschberg());
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::FreeBayesVCFImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::SequenceManipulation::SequenceManipulation() : sequence_manip_impl_ptr_(std::make_unique<kgl::SequenceManipulation::SequenceManipImpl>()) {}
kgl::SequenceManipulation::~SequenceManipulation() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


std::string kgl::SequenceManipulation::compareSequencesDNA(const std::string& reference_str,
                                                        const std::string& compare_str,
                                                        CompareScore_t& score) const {

  return sequence_manip_impl_ptr_->myersHirschberg(reference_str, compare_str, score);

}


std::string kgl::SequenceManipulation::compareSequencesAmino(const std::string& reference_str,
                                                                const std::string& compare_str,
                                                                CompareScore_t& score) const {

  return sequence_manip_impl_ptr_->myersHirschberg(reference_str, compare_str, score);

}


std::string kgl::SequenceManipulation::compareSequencesMultiple(const std::vector<std::string>& sequence_vector) const {

  return sequence_manip_impl_ptr_->multipleAlign(sequence_vector);

}

kgl::CompareScore_t kgl::SequenceManipulation::compareMyerHirschberg(const std::string& reference_str,
                                                                     const std::string& compare_str) const {

  CompareScore_t score;
  sequence_manip_impl_ptr_->myersHirschberg(reference_str, compare_str, score);
  return score;

}
