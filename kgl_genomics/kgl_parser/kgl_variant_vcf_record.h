//
// Created by kellerberrin on 15/4/20.
//

#ifndef KGL_VARIANT_FILE_H
#define KGL_VARIANT_FILE_H

#include "kgl_genome_types.h"
#include <string>
#include <vector>



namespace kellerberrin::genome {   //  organization::project level namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic VCF Record LIne.
////////////////////////////////////////////////////////////////////////////////////////////////////////

class VcfRecord {

public:

  explicit VcfRecord(IOLineRecord&& line_record) noexcept {

    auto [line_count, line_str] = line_record.getLineData();
    line_count_ = line_count;
    line_record_ptr = std::make_unique<std::string>(std::move(line_str));

  }

  std::unique_ptr<const std::string> line_record_ptr;
  size_t line_count_{0};
  bool EOF_{false};
  //
  // Parsed VCF fields
  //
  // Id of the reference sequence.
  ContigId_t contig_id;
  // Position on the reference.
  ContigOffset_t offset{INVALID_POS};
  // Textual identifier of the variant.
  std::string id;
  // Bases in the reference.
  std::string ref;
  // Bases in the alternatives, COMMA-separated.
  std::string alt;
  // Quality
  double qual{MISSING_QUAL};
  // Value of FILTER field.
  std::string filter;
  // Value of INFO field.
  std::string info;
  // Value of FORMAT field.
  std::string format;
  // The genotype infos.
  std::vector<std::string> genotypeInfos;

  // Constant for invalid position.
  static constexpr const ContigOffset_t INVALID_POS = std::numeric_limits<ContigOffset_t>::max();
  // Undefined quality number.
  static constexpr const double MISSING_QUAL = std::numeric_limits<double>::lowest();

  [[nodiscard]] bool EOFRecord() const { return EOF_; }
  // ReturnType the EOF marker.
  [[nodiscard]] static VcfRecord createEOFMarker() { return {}; }

private:

  // Default constructor.
  VcfRecord() { EOF_ = true; }

};


} // namespace

#endif //KGL_VARIANT_FILE_H
