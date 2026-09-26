//
// Created by kellerberrin on 28/09/22.
//

#ifndef KGL_FASTA_H
#define KGL_FASTA_H


#include "kgl_genome_genome.h"
#include "kgl_sequence_view.h"

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>


namespace kellerberrin::genome {   //  organization::project level namespace

// Creates an instance of a Genome database object.
// Parses the input gff(3) file and annotates it with a fasta sequence

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passed in a vector to write a fasta sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////

class WriteFastaSequence {

public:

  WriteFastaSequence() = delete;
  WriteFastaSequence( std::string fasta_id,
                      std::string fasta_description,
                      SequenceRef sequence) : fasta_id_(std::move(fasta_id)),
                                              fasta_description_(std::move(fasta_description)),
                                              fasta_sequence_(std::move(sequence)) {}

  [[nodiscard]] const std::string& fastaId() const { return fasta_id_; }
  [[nodiscard]] const std::string& fastaDescription() const { return fasta_description_; }
  [[nodiscard]] const SequenceRef& fastaSequence() const { return fasta_sequence_; }

private:

  std::string fasta_id_;
  std::string fasta_description_;
  SequenceRef fasta_sequence_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Passed in a vector from fasta sequence read.
//
//////////////////////////////////////////////////////////////////////////////////////////////////

class ReadFastaSequence {

public:

  ReadFastaSequence() = delete;
  ReadFastaSequence(ReadFastaSequence&& moved) = default;
  ReadFastaSequence( std::string&& fasta_id,
                     std::string&& fasta_description,
                     std::string&& fasta_sequence) : fasta_id_(std::move(fasta_id)),
                                                     fasta_description_(std::move(fasta_description)),
                                                     fasta_sequence_ptr_(std::make_unique<std::string>(std::move(fasta_sequence))) {}
  ReadFastaSequence( std::string&& fasta_id,
                     std::string&& fasta_description,
                     std::unique_ptr<std::string>&& fasta_sequence_ptr) : fasta_id_(std::move(fasta_id)),
                                                                          fasta_description_(std::move(fasta_description)),
                                                                          fasta_sequence_ptr_(std::move(fasta_sequence_ptr)) {

    // Defensive; the accessor dereferences the pointer unchecked.
    if (not fasta_sequence_ptr_) {

      fasta_sequence_ptr_ = std::make_unique<std::string>();

    }

  }

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
// Read and write Fasta files.
// Static object to provide data hiding and namespace.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseFasta {

public:

  ParseFasta() = delete;
  ~ParseFasta() = delete;

  [[nodiscard]] static std::shared_ptr<GenomeReference> readFastaFile(const std::string& organism, const std::string& fasta_file_name);

  [[nodiscard]] static bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);

  [[nodiscard]] static bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);


private:

  static constexpr char FASTA_COMMENT_{';'};
  static constexpr char FASTA_ID_{'>'};
  static constexpr size_t FASTA_LINE_LENGTH_{60};

  // Line predicates shared by the reader state machine.
  [[nodiscard]] static bool isSkippableLine(std::string_view line);
  [[nodiscard]] static bool isFastaIdLine(std::string_view line);

  [[nodiscard]] static ReadFastaSequence createFastaSequence( const std::string& fasta_id,
                                                              const std::string& fasta_comment,
                                                              const std::vector<std::string>& fasta_lines);


};



}   // end namespace


#endif //KGL_FASTA_H
