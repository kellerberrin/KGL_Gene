//
// Created by kellerberrin on 28/09/22.
//

#ifndef KGL_GFF3_H
#define KGL_GFF3_H


#include <memory>
#include <string>
#include <map>
#include "kgl_genome_prelim.h"
#include "kgl_genome_genome.h"
#include "kgl_sequence_virtual.h"
#include "kel_logging.h"


namespace kellerberrin::genome {   //  organization::project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Important - the offsets of the features use the half open interval [begin, end) convention and are ZERO based.
// This is different from the gff3 format which are 1-based and are specified as the closed interval [begin, end].
// Thus begin_offset_ = gff.begin - 1
// And end_offset_ = gff.end
// Note that always begin_offset_ < end_offset_. They are not strand adjusted.
// They always correspond to zero based offsets on the relevant contig.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GffRecord {

public:

  GffRecord() = default;
  ~GffRecord() = default;

  [[nodiscard]] bool contig(const std::string_view &contig);
  [[nodiscard]] bool source(const std::string_view &source);
  [[nodiscard]] bool type(const std::string_view &type);
  [[nodiscard]] bool convertStartOffset(const std::string_view &offset);  // This decrements (--offset) to a ZERO based offset.
  [[nodiscard]] bool convertEndOffset(const std::string_view &offset);   // No value change. All offsets are represented as [begin, end).
  [[nodiscard]] bool score(const std::string_view &offset);
  [[nodiscard]] bool phase(const std::string_view &offset);
  [[nodiscard]] bool strand(const std::string_view &strand);
  [[nodiscard]] bool attributes(const std::vector<std::pair<std::string_view, std::string_view>> &tag_value_pairs);

  [[nodiscard]] const std::string &contig() const { return contig_; }
  [[nodiscard]] const std::string &source() const { return source_; }
  [[nodiscard]] const std::string &type() const { return type_; }
  [[nodiscard]] ContigOffset_t begin() const { return begin_position_; }
  [[nodiscard]] ContigOffset_t end() const { return end_position_; }
  [[nodiscard]] double score() const { return score_; }
  [[nodiscard]] uint32_t phase() const { return phase_; }
  [[nodiscard]] StrandSense strand() const { return strand_; }
  [[nodiscard]] const Attributes &attributes() const { return record_attributes_; }

  static constexpr const char *MISSING_VALUE = ".";
  static constexpr ContigOffset_t const INVALID_OFFSET = std::numeric_limits<ContigOffset_t>::max();
  static constexpr double const NO_SCORE = std::numeric_limits<double>::max();
  static constexpr uint32_t const NO_PHASE = 0;
  static constexpr uint32_t const MAX_PHASE = 2;
  static constexpr uint32_t const INVALID_PHASE = std::numeric_limits<uint32_t>::max();

private:

  std::string contig_;
  std::string source_;
  std::string type_;
  ContigOffset_t begin_position_{INVALID_OFFSET};
  ContigOffset_t end_position_{INVALID_OFFSET};
  double score_{NO_SCORE};
  uint32_t phase_{INVALID_PHASE};
  StrandSense strand_{StrandSense::FORWARD};
  Attributes record_attributes_;

  static constexpr const char *STRAND_FORWARD_CHAR{"+"};
  static constexpr const char *STRAND_REVERSE_CHAR{"-"};
  static constexpr std::errc const ERRC_SUCCESS = std::errc();


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read and write Gff3 files.
// Static object to provide data hiding and namespace.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseGff3 {

public:

  ParseGff3() = delete;

  ~ParseGff3() = delete;

  static void readGffFile(const std::string &gff_file_name, GenomeReference &genome_db);

  static std::pair<bool, std::vector<std::unique_ptr<GffRecord>>> readGffFile(const std::string &file_name);

  static std::pair<bool, std::unique_ptr<GffRecord>> parseGff3Record(std::unique_ptr<std::string> &&gff_line);

  static bool parseGffRecord(GenomeReference& genome_db, const GffRecord& gff_record);


private:

  static constexpr const char GFF_COMMENT_{'#'};
  static constexpr const char GFF3_FIELD_DELIM_{'\t'};
  static constexpr const char GFF3_TAG_FIELD_DELIMITER_{';'};
  static constexpr const char GFF3_TAG_ITEM_DELIMITER_{'='};
  static constexpr const size_t GFF3_ITEM_TAG_NAME_{2};

  static constexpr const size_t GFF3_FIELD_COUNT_{9};
  static constexpr const size_t GFF3_CONTIG_FIELD_IDX_{0};
  static constexpr const size_t GFF3_SOURCE_FIELD_IDX_{1};
  static constexpr const size_t GFF3_TYPE_FIELD_IDX_{2};
  static constexpr const size_t GFF3_START_FIELD_IDX_{3};
  static constexpr const size_t GFF3_END_FIELD_IDX_{4};
  static constexpr const size_t GFF3_SCORE_FIELD_IDX_{5};
  static constexpr const size_t GFF3_STRAND_FIELD_IDX_{6};
  static constexpr const size_t GFF3_PHASE_FIELD_IDX_{7};
  static constexpr const size_t GFF3_TAG_FIELD_IDX_{8}; // The "xxxx=yyyy;aaaa=bbbb;..." tag field offset.



};




} // end namespace



#endif //KGL_GFF3_H
