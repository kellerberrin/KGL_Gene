//
// Created by kellerberrin on 22/02/18.
//


#include <iostream>

#include <tuple>                        // for std::make_pair
#include <iterator>
//#include <span>

#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/std/ranges>                    // include all of the standard library's views
#include <seqan3/range/view/all.hpp>            // include all of SeqAn's views
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>

#include <edlib.h>

#include "kgl_sequence_distance_impl.h"
#include "kel_exec_env.h"
#include "kgl_genome_types.h"

namespace kgl = kellerberrin::genome;


class kgl::SequenceDistanceImpl::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;

  kgl::CompareDistance_t LevenshteinGlobalSeqan3(const std::string& sequenceA, const std::string& sequenceB) const;

  kgl::CompareDistance_t LevenshteinGlobal(const std::string& sequenceA, const std::string& sequenceB) const;
  kgl::CompareDistance_t LevenshteinLocal(const std::string& sequenceA, const std::string& sequenceB) const;
  kgl::CompareDistance_t globalblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;
  kgl::CompareDistance_t localblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;

  
};



kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::globalblosum80Distance(const std::string& sequenceA,
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




kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::localblosum80Distance(const std::string& sequenceA,
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

  long score = seqan::localAlignment(align, seqan::Blosum80(open_gap, extend_gap));

  return static_cast<double>(score) * -1.0;  // Invert the scores.

}



kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinGlobalSeqan3(const std::string& sequenceA,
                                                                                             const std::string& sequenceB) const {

  std::vector<seqan3::dna5> Asequence{sequenceA | seqan3::view::char_to<seqan3::dna5>};
  std::vector<seqan3::dna5> Bsequence{sequenceB | seqan3::view::char_to<seqan3::dna5>};

  auto config = seqan3::align_cfg::edit;

  kgl::CompareDistance_t edit_distance = 0.0;

  for (auto const & res : seqan3::align_pairwise(std::tie(Asequence, Bsequence), config))
  {
    edit_distance = std::fabs(res.score());
  }

  return edit_distance;

}


kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinGlobal(const std::string& sequenceA,
                                                                                     const std::string& sequenceB) const {

  EdlibAlignResult result = edlibAlign(sequenceA.c_str(),
                                       sequenceA.size(),
                                       sequenceB.c_str(),
                                       sequenceB.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_LOC, NULL, 0));

  if (result.status != EDLIB_STATUS_OK) {

    ExecEnv::log().error("Problem calculating Global Levenshtein distance using edlib; sequenceA: {}, sequenceB: {}",
                              sequenceA, sequenceB);
    edlibFreeAlignResult(result);
    return 0;
  }

  kgl::CompareDistance_t distance = std::fabs(result.editDistance);
  edlibFreeAlignResult(result);

  return distance;

}


kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinLocal(const std::string& sequenceA,
                                                                                    const std::string& sequenceB) const {



  EdlibAlignResult result;
  // The smaller sequence is presented first with local sequence matching.
  // A distance metric must always be symmetric -> d(x,y) = d(y,x).
  // This could be a bug in EDLIB.
  if (sequenceA.size() <= sequenceB.size()) {

    result = edlibAlign(sequenceA.c_str(),
                        sequenceA.size(),
                        sequenceB.c_str(),
                        sequenceB.size(),
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));


  } else {

    result = edlibAlign(sequenceB.c_str(),
                        sequenceB.size(),
                        sequenceA.c_str(),
                        sequenceA.size(),
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));

  }


  if (result.status != EDLIB_STATUS_OK) {

    ExecEnv::log().error("Problem calculating Local Levenshtein distance using edlib; sequenceA: {}, sequenceB: {}",
                              sequenceA, sequenceB);
    edlibFreeAlignResult(result);
    return 0;
  }

  kgl::CompareDistance_t distance = std::fabs(result.editDistance);

  edlibFreeAlignResult(result);

  kgl::CompareDistance_t compare = LevenshteinGlobalSeqan3(sequenceA, sequenceB);

  ExecEnv::log().info("Local Levenshtein distance using edlib: {}, Global seqan3: {}", distance, compare);

  return distance;


}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto VcfFactory::FreeBayesVCFImpl.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::SequenceDistanceImpl::SequenceDistanceImpl() : sequence_manip_impl_ptr_(std::make_unique<kgl::SequenceDistanceImpl::SequenceManipImpl>()) {}
kgl::SequenceDistanceImpl::~SequenceDistanceImpl() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete pimpl type.


// Distance

kgl::CompareDistance_t kgl::SequenceDistanceImpl::LevenshteinGlobal(const std::string& sequenceA,
                                                                  const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->LevenshteinGlobal(sequenceA, sequenceB);

}


kgl::CompareDistance_t kgl::SequenceDistanceImpl::LevenshteinLocal(const std::string& sequenceA,
                                                                 const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->LevenshteinLocal(sequenceA, sequenceB);

}


kgl::CompareDistance_t kgl::SequenceDistanceImpl::globalblosum80Distance(const std::string& sequenceA,
                                                                       const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->globalblosum80Distance(sequenceA, sequenceB);

}



kgl::CompareDistance_t kgl::SequenceDistanceImpl::localblosum80Distance(const std::string& sequenceA,
                                                                         const std::string& sequenceB) const {

  return sequence_manip_impl_ptr_->localblosum80Distance(sequenceA, sequenceB);

}

