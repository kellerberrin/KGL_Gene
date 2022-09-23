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
  static constexpr const size_t GFF3_FIELD_COUNT_{9};

  static void readGffFile( const std::string &gff_file_name, GenomeReference& genome_db);
  static void readGffFile(const std::string& file_name);
  static ReadFastaSequence createFastaSequence( const std::string& fasta_id,
                                                const std::string& fasta_comment,
                                                const std::vector<std::unique_ptr<std::string>>& fasta_lines);

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Important - the offsets of the features use the half open interval [begin, end) convention and are ZERO based.
// This is different from the gff format which are 1-based and are specified as the closed interval [begin, end].
// Thus begin_offset_ = gff.begin - 1
// And end_offset_ = gff.end
// Note that always begin_offset_ < end_offset_. They are not strand adjusted.
// They always correspond to zero based offsets on the relevant contig.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GffRecord
{

  static constexpr ContigOffset_t const INVALID_POS = std::numeric_limits<ContigOffset_t>::max();

  std::string ref_;
  std::string source_;
  std::string type_;
  std::map<std::string, std::string> tagNamesValues_;
  ContigOffset_t beginPos_{INVALID_POS};
  ContigOffset_t endPos_{INVALID_POS};
  FeatureSequence feature_sequence_;

};





}   // end namespace


#endif