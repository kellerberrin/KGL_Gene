//
// Created by kellerberrin on 22/02/18.
//


// #define SEQAN3 1
#ifdef SEQAN3

#include <seqan3/alphabet/all.hpp>
#include <seqan3/alignment/all.hpp>
#include <seqan3/alphabet/all.hpp>

#endif

#include "edlib.h"

#include "kgl_sequence_distance_impl.h"
#include "kel_exec_env.h"

#include <cmath>
#include <cstring>

namespace kgl = kellerberrin::genome;


kgl::CompareDistance_t kgl::LevenshteinGlobalImpl(const char* sequenceA,
                                                  size_t sequenceA_size,
                                                  const char* sequenceB,
                                                  size_t sequenceB_size) {

  EdlibAlignResult result = edlibAlign(sequenceA,
                                       sequenceA_size,
                                       sequenceB,
                                       sequenceB_size,
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



kgl::CompareDistance_t kgl::LevenshteinLocalImpl(const char* sequenceA,
                                                 size_t sequenceA_size,
                                                 const char* sequenceB,
                                                 size_t sequenceB_size) {


  EdlibAlignResult result;
  // The smaller sequence is presented first with local sequence matching.
  // A calculateDistance metric must always be symmetric -> d(x,y) = d(y,x).
  // This could be a bug in EDLIB.

  if (sequenceA_size <= sequenceB_size) {

    result = edlibAlign(sequenceA,
                        sequenceA_size,
                        sequenceB,
                        sequenceB_size,
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));


  } else {

    result = edlibAlign(sequenceB,
                        sequenceB_size,
                        sequenceA,
                        sequenceA_size,
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


