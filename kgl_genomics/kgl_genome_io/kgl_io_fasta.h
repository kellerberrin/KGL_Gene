//
// Created by kellerberrin on 28/09/22.
//

#ifndef KGL_FASTA_H
#define KGL_FASTA_H


#include "kgl_genome_prelim.h"
#include "kgl_genome_genome.h"
#include "kgl_sequence_virtual.h"

#include <memory>
#include <string>
#include <map>


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
// Read and write Fasta files.
// Static object to provide data hiding and namespace.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseFasta {

public:

  ParseFasta() =delete;
  ~ParseFasta() = delete;

  [[nodiscard]] static std::shared_ptr<GenomeReference> readFastaFile(const std::string& organism, const std::string& fasta_file_name);

  [[nodiscard]] static bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);

  [[nodiscard]] static bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);


private:

  static constexpr const char FASTA_COMMENT_{';'};
  static constexpr const char FASTA_ID_{'>'};

  static ReadFastaSequence createFastaSequence( const std::string& fasta_id,
                                                const std::string& fasta_comment,
                                                const std::vector<std::unique_ptr<const std::string>>& fasta_lines);


};



}   // end namespace


#endif //KGL_FASTA_H
