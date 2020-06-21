//
// Created by kellerberrin on 5/5/20.
//

#ifndef KGL_ANALYSIS_ALL_H
#define KGL_ANALYSIS_ALL_H

// Includes all of the analytic header files for convenience.

#include "kgl_analysis_null.h" // The base class
#include "kgl_analysis_interval.h" // basic sequence and variant statistics.
#include "kgl_analysis_info_filter.h" // basic sequence and variant statistics.

namespace kgl = kellerberrin::genome;

extern "C" {


kgl::FactoryAnalysisVector kgl_packageAnalysisFactory();


}

#endif //KGL_KGL_ANALYSIS_ALL_H
