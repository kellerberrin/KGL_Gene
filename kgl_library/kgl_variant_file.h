//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H

#include "kgl_genome_types.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <string>
#include <vector>
#include <fstream>

namespace bio = boost::iostreams;
namespace kellerberrin::genome {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic VCF Record LIne.
////////////////////////////////////////////////////////////////////////////////////////////////////////

class VcfRecord
{
public:

  // Default constructor.
  VcfRecord() : offset(INVALID_POS), qual(MISSING_QUAL) {}

  static std::unique_ptr<VcfRecord> EOF_RECORD() { return std::unique_ptr<VcfRecord>(); }

  // Numeric id of the reference sequence.
  ContigId_t contig_id;
  // Position on the reference.
  ContigOffset_t offset;
  // Textual identifier of the variant.
  std::string id;
  // Bases in the reference.
  std::string ref;
  // Bases in the alternatives, comma-separated.
  std::string alt;
  // Quality
  double qual;
  // Value of FILTER field.
  std::string filter;
  // Value of INFO field.
  std::string info;
  // Value of FORMAT field.
  std::string format;
  // The genotype infos.
  std::vector<std::string> genotypeInfos;


private:
  // Constant for invalid position.
  static constexpr const ContigOffset_t INVALID_POS = std::numeric_limits<ContigOffset_t>::max();
  // Undefined quality number.
  static constexpr const double MISSING_QUAL = std::numeric_limits<double>::lowest();

};


} // end namespace


#endif //KGL_VARIANT_FILE_H
