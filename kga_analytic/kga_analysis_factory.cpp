//
// Created by kellerberrin on 5/5/20.
//


#include "kga_analysis_null.h" // The "do nothing" class and a template for additional analysis objects.
#include "kga_analysis_sequence.h" // Does "correctness" verification on any specified data files (duplicate variants etc).
#include "kga_analysis_json.h" // Processes Allele Json files and writes resultant the information to file.
#include "kga_analysis_interval.h" // Basic sequence and variant statistics.
#include "kga_analysis_info_filter.h" // Age related variant statistics, only applicable to Gnomad variant data.
#include "kga_analysis_inbreed.h" // Analysis of Diploid phased population read from the 1000 genomes project.
#include "kga_analysis_PfEMP.h"   // Analyze the P. Falciparum protein family for different Pf Genomes.
#include "kga_analysis_mutation.h"   // Analyze mutations in genes and genomic regions.
#include "kga_analysis_literature.h"   // Analyze Pubmed literature.


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Lookup table to dynamically create analysis objects specified in package XML files.
// A list of available analytic objects passed into the PackageAnalysis object.
// Any analysis specified in the Package XML must match the std::string IDENT
// of one of these objects.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kgl = kellerberrin::genome;
namespace kga = kellerberrin::genome::analysis;


kgl::VirtualAnalysis::AnalysisFactoryMap const kgl::VirtualAnalysis::analysis_factory_map_ = {

{ kga::NullAnalysis::IDENT, kga::NullAnalysis::factory },
{ kga::SequenceAnalysis::IDENT, kga::SequenceAnalysis::factory },
{ kga::IntervalAnalysis::IDENT, kga::IntervalAnalysis::factory },
{ kga::InfoFilterAnalysis::IDENT, kga::InfoFilterAnalysis::factory },
{ kga::InbreedAnalysis::IDENT, kga::InbreedAnalysis::factory },
{ kga::PfEMPAnalysis::IDENT, kga::PfEMPAnalysis::factory },
{ kga::MutationAnalysis::IDENT, kga::MutationAnalysis::factory },
{ kga::JsonAnalysis::IDENT, kga::JsonAnalysis::factory },
{ kga::LiteratureAnalysis::IDENT, kga::LiteratureAnalysis::factory }

};


