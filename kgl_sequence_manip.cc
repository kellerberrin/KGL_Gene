//
// Created by kellerberrin on 4/01/18.
//

#include <iostream>
#include <seqan/align.h>

#include "kgl_sequence_manip.h"


namespace kgl = kellerberrin::genome;
using namespace seqan;



class kgl::SequenceManipulation::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;

  std::string compareSequences(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

  std::string gap_example(const std::string& reference_str, const std::string& compare_str) const;

  std::string graph_example(const std::string& reference_str, const std::string& compare_str) const;

  std::string align_example(const std::string& reference_str, const std::string& compare_str, CompareScore_t& score) const;

private:


};



std::string kgl::SequenceManipulation::SequenceManipImpl::compareSequences(const std::string& reference_str,
                                                                           const std::string& compare_str,
                                                                           CompareScore_t& score) const {

  return align_example(reference_str, compare_str, score);

}


std::string kgl::SequenceManipulation::SequenceManipImpl::align_example(const std::string& reference_str,
                                                                        const std::string& compare_str,
                                                                        CompareScore_t& score) const {
  typedef String<char> TSequence;
  typedef Align<TSequence, ArrayGaps> TAlign;
//  typedef Row<TAlign>::Type TRow;
//  typedef Iterator<TRow>::Type TRowIterator;

//  TSequence seq1 = "AAGUGACUUAUUG";
//  TSequence seq2 = "AGUCGGAUCUACUG";

  TSequence seq1 = reference_str.c_str();
  TSequence seq2 = compare_str.c_str();

  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq2);

  std::stringstream ss;

  score = globalAlignment(align, MyersHirschberg());
  ss << align << std::endl;

  return ss.str();

}



std::string kgl::SequenceManipulation::SequenceManipImpl::graph_example(const std::string& reference_str,
                                                                           const std::string&) const {

  typedef String<char> TSequence;
  typedef StringSet<TSequence> TStringSet;
  typedef StringSet<TSequence, Dependent<> > TDepStringSet;
  typedef Graph<Alignment<TDepStringSet> > TAlignGraph;

  // Initializing the sequences and the string set.
  TSequence seq1 = "GARFIELDTHECAT";
  TSequence seq2 = "GARFIELDTHEBIGCAT";
  TSequence seq3 = "THEBIGCAT";

  TStringSet strings;
  appendValue(strings, seq1);
  appendValue(strings, seq2);
  appendValue(strings, seq3);

  // Load the string set into the Alignment Graph.
  TAlignGraph alignG(strings);

  // Add two vertices covering "GARFIELD" in the first and the second sequence and connect them with an edge.
  addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG), 0), 0, 8),
          addVertex(alignG, positionToId(stringSet(alignG), 1), 0, 8));

  // Add two vertices covering "THE" in the first and the second sequence and connect them with an edge.
  addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG), 0), 8, 3),
          addVertex(alignG, positionToId(stringSet(alignG), 1), 8, 3));

  // Find the vertex covering "THE" in the first sequence and add the vertex covering "THE" in the third sequence and connect them with an edge.
  addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG), 0), 8),
          addVertex(alignG, positionToId(stringSet(alignG), 2), 0, 3));

  // Find the vertices covering "THE" in the second and the third sequence and connect them with an edge.
  addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG), 1), 8),
          findVertex(alignG, positionToId(stringSet(alignG), 2), 0));

  // Add two vertices covering "FAT" in the second and the third sequence and connect them with an edge.
  addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG), 1), 11, 3),
          addVertex(alignG, positionToId(stringSet(alignG), 2), 3, 3));

  // Add two vertices covering "CAT" in the first and the second sequence and connect them with an edge.
  addEdge(alignG, addVertex(alignG, positionToId(stringSet(alignG), 0), 11, 3),
          addVertex(alignG, positionToId(stringSet(alignG), 1), 14, 3));

  // Find the vertex covering "CAT" in the first sequence and add the vertex covering "CAT" in the third sequence and connect them with an edge.
  addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG), 0), 11),
          addVertex(alignG, positionToId(stringSet(alignG), 2), 6, 3));

  // Find the vertices covering "CAT" in the second and the third sequence and connect them with an edge.
  addEdge(alignG, findVertex(alignG, positionToId(stringSet(alignG), 1), 14),
          findVertex(alignG, positionToId(stringSet(alignG), 2), 6));

  std::cout << alignG << std::endl;

  return reference_str;

}

std::string kgl::SequenceManipulation::SequenceManipImpl::gap_example(const std::string& reference_str,
                                                                           const std::string&) const {

  // Defining all types that are needed.
  using TSequence = seqan::String<char>;
  using TAlign = seqan::Align<TSequence, seqan::ArrayGaps>;
  using TRow = seqan::Row<TAlign>::Type;
  using TRowIterator = seqan::Iterator<TRow>::Type ;

  TSequence seq1 = "ACGTCACCTC";
  TSequence seq2 = "ACGGGCCTATC";

//  TSequence seq1 = reference_str.c_str();
//  TSequence seq2 = compare_str.c_str();


  // Initializing the align object.
  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq2);

  std::cout << " align before gap:" << align << std::endl;
  // Use references to the rows of align.
  TRow & row1 = row(align, 0);
  TRow & row2 = row(align, 1);

  // Insert gaps.
  insertGaps(row1, 2, 2);
  insertGap(row1, 7);  // We need to pass the view position which is changed due to the previous insertion.
  insertGaps(row2, 9, 2);
  std::cout << " align after gap:" << align << std::endl;

  std::cout << std::endl << "ViewToSource1: ";
  for (unsigned i = 0; i < length(row1); ++i)
    std::cout << toSourcePosition(row1, i) << ",";

  std::cout << std::endl << "ViewToSource2: ";
  for (unsigned i = 0; i < seqan::length(row2); ++i)
    std::cout << seqan::toSourcePosition(row2, i) << ",";
  std::cout << std::endl;

  std::cout << std::endl << "SourceToView1: ";
  for (unsigned i = 0; i < seqan::length(seqan::source(row1)); ++i)
    std::cout << seqan::toViewPosition(row1, i) << ",";

  std::cout << std::endl << "SourceToView2: ";
  for (unsigned i = 0; i < seqan::length(seqan::source(row2)); ++i)
    std::cout << seqan::toViewPosition(row2, i) << ",";
  std::cout << std::endl;


  // Initialize the row iterators.
  TRowIterator itRow1 = seqan::begin(row1);
  TRowIterator itEndRow1 = seqan::end(row1);
  TRowIterator itRow2 = seqan::begin(row2);

  // Iterate over both rows simultaneously.
  int gapCount = 0;
  for (; itRow1 != itEndRow1; ++itRow1, ++itRow2)
  {
    if (seqan::isGap(itRow1))
    {
      gapCount += seqan::countGaps(itRow1);  // Count the number of consecutive gaps from the current position in row1.
      itRow1 += seqan::countGaps(itRow1);    // Jump to next position to check for gaps.
      itRow2 += seqan::countGaps(itRow1);    // Jump to next position to check for gaps.
    }
    if (seqan::isGap(itRow2))
    {
      gapCount += seqan::countGaps(itRow2);  // Count the number of consecutive gaps from the current position in row2.
      itRow1 += seqan::countGaps(itRow2);    // Jump to next position to check for gaps.
      itRow2 += seqan::countGaps(itRow2);    // Jump to next position to check for gaps.
    }
  }
  // Print the result.
  std::cout << "Number of gaps: " << gapCount << std::endl;

  return reference_str;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::VcfFileImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::SequenceManipulation::SequenceManipulation() : sequence_manip_impl_ptr_(std::make_unique<kgl::SequenceManipulation::SequenceManipImpl>()) {}
kgl::SequenceManipulation::~SequenceManipulation() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


std::string kgl::SequenceManipulation::compareSequences(const std::string& reference_str,
                                                        const std::string& compare_str,
                                                        CompareScore_t& score) const {

  return sequence_manip_impl_ptr_->compareSequences(reference_str, compare_str, score);

}

