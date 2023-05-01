//
// Created by kellerberrin on 22/02/18.
//


#include <iostream>

#include <tuple>                        // for std::make_pair
#include <iterator>
#include <string>

// #define SEQAN3 1
#ifdef SEQAN3

#include <seqan3/alphabet/all.hpp>
#include <seqan3/alignment/all.hpp>
#include <seqan3/alphabet/all.hpp>

#endif

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




#ifdef SEQAN3

kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::globalblosum80Distance(const std::string& sequenceA,
                                                                                          const std::string& sequenceB) const {

  using TSequence = seqan3::String<seqan3::AminoAcid> ;
  using TAlign = seqan3::Align<TSequence, seqan3::ArrayGaps> ;
  int open_gap = -8;
  int extend_gap = -3;

  TSequence seq1 = sequenceA.c_str();
  TSequence seq2 = sequenceB.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan3::assignSource(row(align, 0), seq1);
  seqan3::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  long score = seqan3::globalAlignment(align, seqan3::Blosum80(open_gap, extend_gap));

  return static_cast<double>(score) * -1.0;  // Invert the scores.

}

#else

kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::globalblosum80Distance(const std::string&,
                                                                                            const std::string&) const {
  return 0.0;

}

#endif





#ifdef SEQAN3

kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::localblosum80Distance(const std::string& sequenceA,
                                                                                            const std::string& sequenceB) const {

  using TSequence = seqan3::String<seqan3::AminoAcid> ;
  using TAlign = seqan3::Align<TSequence, seqan3::ArrayGaps> ;
  int open_gap = -8;
  int extend_gap = -3;

  TSequence seq1 = sequenceA.c_str();
  TSequence seq2 = sequenceB.c_str();

  TAlign align;
  resize(rows(align), 2);
  seqan3::assignSource(row(align, 0), seq1);
  seqan3::assignSource(row(align, 1), seq2);

  std::stringstream ss;

  long score = seqan3::localAlignment(align, seqan3::Blosum80(open_gap, extend_gap));

  return static_cast<double>(score) * -1.0;  // Invert the scores.

}

#else

kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::localblosum80Distance(const std::string&,
                                                                                           const std::string&) const {

  return 0.0;

}

#endif



#ifdef SEQAN3

  kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinGlobalSeqan3(const std::string& sequenceA,
                                                                                             const std::string& sequenceB) const {

  auto Asequence = seqan3::views::char_to<seqan3::dna5>(sequenceA);
  auto Bsequence = seqan3::views::char_to<seqan3::dna5>(sequenceB);

  auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

  kgl::CompareDistance_t edit_distance = 0.0;
  for (auto const & res : seqan3::align_pairwise(std::tie(Asequence, Bsequence), config)) {

    edit_distance = std::fabs(res.score());

  }

  return edit_distance;

}

#else

  kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinGlobalSeqan3(const std::string&,
                                                                                               const std::string&) const {


    return 0.0;

  }

#endif


kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinGlobal(const std::string& sequenceA,
                                                                                     const std::string& sequenceB) const {

  EdlibAlignResult result = edlibAlign(sequenceA.c_str(),
                                       sequenceA.size(),
                                       sequenceB.c_str(),
                                       sequenceB.size(),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_LOC, NULL, 0));

  if (result.status != EDLIB_STATUS_OK) {

    ExecEnv::log().error("Problem calculating Global Levenshtein calculateDistance using edlib; sequenceA: {}, sequenceB: {}",
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
  // A calculateDistance metric must always be symmetric -> d(x,y) = d(y,x).
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

    ExecEnv::log().error("Problem calculating Local Levenshtein calculateDistance using edlib; sequenceA: {}, sequenceB: {}",
                              sequenceA, sequenceB);
    edlibFreeAlignResult(result);
    return 0;
  }

  kgl::CompareDistance_t distance = std::fabs(result.editDistance);

  edlibFreeAlignResult(result);

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

