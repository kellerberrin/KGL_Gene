//
// kgl_sequence_distance_impl.cpp — edlib backed Levenshtein distance.
//
// Created by kellerberrin on 22/02/18.
//


#include "edlib.h"

#include "kgl_sequence_distance_impl.h"
#include "kel_exec_env.h"

#include <cmath>


namespace kellerberrin::genome {   //  organization level namespace


CompareDistance_t LevenshteinGlobalImpl(const char* sequenceA,
                                        size_t sequenceA_size,
                                        const char* sequenceB,
                                        size_t sequenceB_size) {

  EdlibAlignResult result = edlibAlign(sequenceA,
                                       static_cast<int32_t>(sequenceA_size),
                                       sequenceB,
                                       static_cast<int32_t>(sequenceB_size),
                                       edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));

  if (result.status != EDLIB_STATUS_OK) {

    ExecEnv::log().error("Problem calculating Global Levenshtein distance using edlib; sequenceA size: {}, sequenceB size: {}",
                         sequenceA_size, sequenceB_size);
    edlibFreeAlignResult(result);
    return 0;

  }

  CompareDistance_t distance = std::fabs(static_cast<double>(result.editDistance));
  edlibFreeAlignResult(result);

  return distance;

}


CompareDistance_t LevenshteinLocalImpl(const char* sequenceA,
                                      size_t sequenceA_size,
                                      const char* sequenceB,
                                      size_t sequenceB_size) {

  // The smaller sequence is presented first with local sequence matching.
  // A distance metric must always be symmetric -> d(x,y) = d(y,x).
  // This could be a bug in EDLIB.

  EdlibAlignResult result;
  if (sequenceA_size <= sequenceB_size) {

    result = edlibAlign(sequenceA,
                        static_cast<int32_t>(sequenceA_size),
                        sequenceB,
                        static_cast<int32_t>(sequenceB_size),
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, nullptr, 0));

  } else {

    result = edlibAlign(sequenceB,
                        static_cast<int32_t>(sequenceB_size),
                        sequenceA,
                        static_cast<int32_t>(sequenceA_size),
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, nullptr, 0));

  }

  if (result.status != EDLIB_STATUS_OK) {

    ExecEnv::log().error("Problem calculating Local Levenshtein distance using edlib; sequenceA size: {}, sequenceB size: {}",
                         sequenceA_size, sequenceB_size);
    edlibFreeAlignResult(result);
    return 0;

  }

  CompareDistance_t distance = std::fabs(static_cast<double>(result.editDistance));

  edlibFreeAlignResult(result);

  return distance;

}


}   // end namespace