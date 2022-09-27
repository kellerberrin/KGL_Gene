///
// Created by kellerberrin on 3/10/17.
//

#ifndef KGL_GFF_FASTA_H
#define KGL_GFF_FASTA_H


#include <memory>
#include <string>
#include <map>
#include "kgl_genome_prelim.h"
#include "kgl_genome_genome.h"
#include "kgl_sequence_virtual.h"
#include "kel_logging.h"


namespace kellerberrin::genome {   //  organization::project level namespace

// Creates an instance of a Genome database object.
// Parses the input gff(3) file and annotates it with a fasta sequence
// This class is a facade; the file reading details are handled by a 3rd party library (Seqan).

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passed in a vector to write a fasta sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////


class WriteFastaSequence {

public:

  WriteFastaSequence() = delete;
  WriteFastaSequence(const WriteFastaSequence& copy) = default;
  WriteFastaSequence( std::string fasta_id,
                      std::string fasta_description,
                      std::shared_ptr<const VirtualSequence> fasta_sequence_ptr) : fasta_id_(std::move(fasta_id)),
                                                                                   fasta_description_(std::move(fasta_description)),
                                                                                   fasta_sequence_ptr_(std::move(fasta_sequence_ptr)) {}
  ~WriteFastaSequence() = default;

  [[nodiscard]] const std::string& fastaId() const { return fasta_id_; }
  [[nodiscard]] const std::string& fastaDescription() const { return fasta_description_; }
  [[nodiscard]] const std::shared_ptr<const VirtualSequence>& fastaSequence() const { return fasta_sequence_ptr_; }

private:

  std::string fasta_id_;
  std::string fasta_description_;
  std::shared_ptr<const VirtualSequence> fasta_sequence_ptr_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passed in a vector from fasta sequence read.
//
//////////////////////////////////////////////////////////////////////////////////////////////////

class ReadFastaSequence {

public:

  ReadFastaSequence() = delete;
  ReadFastaSequence(ReadFastaSequence&& copy) = default;
  ReadFastaSequence( std::string&& fasta_id,
                     std::string&& fasta_description,
                     std::string&& fasta_sequence) : fasta_id_(fasta_id),
                                                     fasta_description_(fasta_description),
                                                     fasta_sequence_ptr_(std::make_unique<std::string>(fasta_sequence)) {}
  ReadFastaSequence( std::string&& fasta_id,
                     std::string&& fasta_description,
                     std::unique_ptr<std::string>&& fasta_sequence_ptr) : fasta_id_(std::move(fasta_id)),
                                                     fasta_description_(std::move(fasta_description)),
                                                     fasta_sequence_ptr_(std::move(fasta_sequence_ptr)) {}

  ~ReadFastaSequence() = default;

  [[nodiscard]] const std::string& fastaId() const { return fasta_id_; }
  [[nodiscard]] const std::string& fastaDescription() const { return fasta_description_; }
  [[nodiscard]] const std::string& fastaSequence() const { return *fasta_sequence_ptr_; }

private:

  std::string fasta_id_;
  std::string fasta_description_;
  std::unique_ptr<std::string> fasta_sequence_ptr_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Important - the offsets of the features use the half open interval [begin, end) convention and are ZERO based.
// This is different from the gff3 format which are 1-based and are specified as the closed interval [begin, end].
// Thus begin_offset_ = gff.begin - 1
// And end_offset_ = gff.end
// Note that always begin_offset_ < end_offset_. They are not strand adjusted.
// They always correspond to zero based offsets on the relevant contig.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GffRecord
{

public:

  GffRecord() = default;
  ~GffRecord() = default;

  [[nodiscard]] bool id(const std::string_view& id);
  [[nodiscard]] bool source(const std::string_view& source);
  [[nodiscard]] bool type(const std::string_view& type);
  [[nodiscard]] bool convertStartOffset(const std::string_view& offset);  // This decrements (--offset) to a ZERO based offset.
  [[nodiscard]] bool convertEndOffset(const std::string_view& offset);   // No value change. All offsets are represented as [begin, end).
  [[nodiscard]] bool score(const std::string_view& offset);
  [[nodiscard]] bool phase(const std::string_view& offset);
  [[nodiscard]] bool strand(const std::string_view& strand);
  [[nodiscard]] bool tagValues(const std::vector<std::pair<std::string_view, std::string_view>>& tag_value_pairs);

  [[nodiscard]] const std::string& id() const { return id_; }
  [[nodiscard]] const std::string& source() const { return source_; }
  [[nodiscard]] const std::string& type() const { return type_; }
  [[nodiscard]] ContigOffset_t begin() const { return begin_position_; }
  [[nodiscard]] ContigOffset_t end() const { return end_position_; }
  [[nodiscard]] double score() const { return score_; }
  [[nodiscard]] uint32_t phase() const { return phase_; }
  [[nodiscard]] StrandSense strand() const { return strand_; }
  [[nodiscard]] const std::map<std::string, std::string>& tagValues() const { return tagValues_; }

  static constexpr const char* MISSING_VALUE = ".";
  static constexpr ContigOffset_t const INVALID_OFFSET = std::numeric_limits<ContigOffset_t>::max();
  static constexpr double const NO_SCORE = std::numeric_limits<double>::max();
  static constexpr uint32_t const NO_PHASE = 0;
  static constexpr uint32_t const MAX_PHASE = 2;
  static constexpr uint32_t const INVALID_PHASE = std::numeric_limits<uint32_t>::max();
  static constexpr const char* STRAND_FORWARD_CHAR{"+"};
  static constexpr const char* STRAND_REVERSE_CHAR{"-"};

  static constexpr std::errc const ERRC_SUCCESS = std::errc();

private:

  std::string id_;
  std::string source_;
  std::string type_;
  ContigOffset_t begin_position_{INVALID_OFFSET};
  ContigOffset_t end_position_{INVALID_OFFSET};
  double score_{NO_SCORE};
  uint32_t phase_{INVALID_PHASE};
  StrandSense strand_{StrandSense::FORWARD};
  std::map<std::string, std::string> tagValues_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Read and write Gff3 abd Fasta files.
// Static object to provide data hiding and namespace.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseGffFasta {

public:

  ParseGffFasta() =delete;
  ~ParseGffFasta() = delete;

  [[nodiscard]] static std::shared_ptr<GenomeReference> readFastaFile(const std::string& organism, const std::string& fasta_file_name);

  [[nodiscard]] static bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);

  [[nodiscard]] static bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);

  [[nodiscard]] static std::shared_ptr<GenomeReference> readFastaGffFile( const std::string& organism,
                                                                          const std::string& fasta_file_name,
                                                                          const std::string& gff_file_name);

private:

  static constexpr const char FASTA_COMMENT_{';'};
  static constexpr const char FASTA_ID_{'>'};
  static constexpr const char GFF_COMMENT_{'#'};
  static constexpr const char GFF3_FIELD_DELIM_{'\t'};
  static constexpr const char GFF3_TAG_FIELD_DELIMITER_{';'};
  static constexpr const char GFF3_TAG_ITEM_DELIMITER_{'='};
  static constexpr const size_t GFF3_ITEM_TAG_NAME_{2};

  static constexpr const size_t GFF3_FIELD_COUNT_{9};
  static constexpr const size_t GFF3_ID_FIELD_IDX_{0};
  static constexpr const size_t GFF3_SOURCE_FIELD_IDX_{1};
  static constexpr const size_t GFF3_TYPE_FIELD_IDX_{2};
  static constexpr const size_t GFF3_START_FIELD_IDX_{3};
  static constexpr const size_t GFF3_END_FIELD_IDX_{4};
  static constexpr const size_t GFF3_SCORE_FIELD_IDX_{5};
  static constexpr const size_t GFF3_STRAND_FIELD_IDX_{6};
  static constexpr const size_t GFF3_PHASE_FIELD_IDX_{7};
  static constexpr const size_t GFF3_TAG_FIELD_IDX_{8}; // The "xxxx=yyyy;aaaa=bbbb;..." tag field offset.

  static void readGffFile( const std::string &gff_file_name, GenomeReference& genome_db);
  static void readGffFile(const std::string& file_name);
  static ReadFastaSequence createFastaSequence( const std::string& fasta_id,
                                                const std::string& fasta_comment,
                                                const std::vector<std::unique_ptr<std::string>>& fasta_lines);
  static std::pair<bool, std::unique_ptr<GffRecord>> parseGff3Record(std::unique_ptr<std::string>&& gff_line);


};





}   // end namespace


#endif