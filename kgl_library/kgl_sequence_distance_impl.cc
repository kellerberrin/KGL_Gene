//
// Created by kellerberrin on 22/02/18.
//


#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <edlib.h>


#include "kgl_sequence_distance_impl.h"
#include "kgl_exec_env.h"
#include "kgl_genome_types.h"


namespace kgl = kellerberrin::genome;


class kgl::SequenceDistanceImpl::SequenceManipImpl {

public:

  SequenceManipImpl() = default;
  ~SequenceManipImpl() = default;
  
  kgl::CompareDistance_t LevenshteinGlobal(const std::string& sequenceA, const std::string& sequenceB) const;
  kgl::CompareDistance_t LevenshteinLocal(const std::string& sequenceA, const std::string& sequenceB) const;
  kgl::CompareDistance_t globalblosum80Distance(const std::string& sequenceA, const std::string& sequenceB) const;

  
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


kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinGlobal(const std::string& sequenceA,
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


kgl::CompareDistance_t kgl::SequenceDistanceImpl::SequenceManipImpl::LevenshteinLocal(const std::string& sequenceA,
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

