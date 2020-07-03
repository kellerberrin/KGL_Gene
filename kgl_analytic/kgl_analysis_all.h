//
// Created by kellerberrin on 5/5/20.
//

#ifndef KGL_ANALYSIS_ALL_H
#define KGL_ANALYSIS_ALL_H

// Includes all of the analytic header files for convenience.

#include "kgl_analysis_virtual.h" // The pure virtual interface definition and base class
#include "kgl_analysis_null.h" // The "do nothing" class and a template for additional analysis objects.
#include "kgl_analysis_interval.h" // Basic sequence and variant statistics.
#include "kgl_analysis_info_filter.h" // Age related variant statistics, only applicable to Gnomad variant data.
#include "kgl_analysis_mutation.h" // Analysis of Diploid phased population read from the 1000 genomes project.


namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A list of available analytic objects passed into the PackageAnalysis object.
// Any analysis specified in the Package XML must match the std::string ident() function
// of one of these objects.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Analytic classes are called virtually and provided with data (reference and variant).
using VirtualAnalysisVector = std::vector<std::unique_ptr<VirtualAnalysis>>;

inline VirtualAnalysisVector getAnalysisVector() {

  VirtualAnalysisVector analysis_vector;

  analysis_vector.push_back(std::make_unique<NullAnalysis>());
  analysis_vector.push_back(std::make_unique<IntervalAnalysis>());
  analysis_vector.push_back(std::make_unique<InfoFilterAnalysis>());
  analysis_vector.push_back(std::make_unique<MutationAnalysis>());

  return analysis_vector;

}




} //namespace



#endif //KGL_KGL_ANALYSIS_ALL_H
