//
// Created by kellerberrin on 5/5/20.
//

#ifndef KGL_ANALYSIS_ALL_H
#define KGL_ANALYSIS_ALL_H

// Includes all of the analytic header files for convenience.

#include "kga_analysis_virtual.h" // The pure virtual interface definition and base class
#include "kga_template_analysis/kga_analysis_null.h" // The "do nothing" class and a template for additional analysis objects.
#include "kga_analysis_sequence.h" // Does "correctness" verification on any specified data files (duplicate variants etc).
#include "kga_analysis_json.h" // Processes Allele Json files and writes resultant the information to file.
#include "kga_analysis_interval.h" // Basic sequence and variant statistics.
#include "kga_analysis_info_filter.h" // Age related variant statistics, only applicable to Gnomad variant data.
#include "kga_analysis_inbreed.h" // Analysis of Diploid phased population read from the 1000 genomes project.
#include "kga_analysis_PfEMP.h"   // Analyze the P. Falciparum protein family for different Pf Genomes.
#include "kga_analysis_mutation.h"   // Analyze mutations in genes and genomic regions.
#include "kga_analysis_literature.h"   // Analyze Pubmed literature.


namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A list of available analytic objects passed into the PackageAnalysis object.
// Any analysis specified in the Package XML must match the std::string ident() function
// of one of these objects.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Analytic classes are called virtually and provided with data and runtime parameters specified in the XML
// definition files.
using VirtualAnalysisVector = std::vector<std::unique_ptr<VirtualAnalysis>>;

inline VirtualAnalysisVector getAnalysisVector() {

  VirtualAnalysisVector analysis_vector;

  analysis_vector.push_back(std::make_unique<NullAnalysis>());
  analysis_vector.push_back(std::make_unique<VerifyAnalysis>());
  analysis_vector.push_back(std::make_unique<IntervalAnalysis>());
  analysis_vector.push_back(std::make_unique<InfoFilterAnalysis>());
  analysis_vector.push_back(std::make_unique<InbreedAnalysis>());
  analysis_vector.push_back(std::make_unique<PfEMPAnalysis>());
  analysis_vector.push_back(std::make_unique<MutationAnalysis>());
  analysis_vector.push_back(std::make_unique<JsonAnalysis>());
  analysis_vector.push_back(std::make_unique<LiteratureAnalysis>());

  return analysis_vector;

}




} //namespace



#endif //KGL_KGL_ANALYSIS_ALL_H
