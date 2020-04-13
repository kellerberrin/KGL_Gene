///
// Created by kellerberrin on 3/10/17.
//

#ifndef KGL_GFF_FASTA_H
#define KGL_GFF_FASTA_H


#include <memory>
#include <string>
#include <map>
#include "kgl_genome_db.h"
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
  ~ReadFastaSequence() = default;

  [[nodiscard]] const std::string& fastaId() const { return fasta_id_; }
  [[nodiscard]] const std::string& fastaDescription() const { return fasta_description_; }
  [[nodiscard]] const std::string& fastaSequence() const { return *fasta_sequence_ptr_; }

private:

  std::string fasta_id_;
  std::string fasta_description_;
  std::unique_ptr<std::string> fasta_sequence_ptr_;

};




class ParseGffFasta {

public:

  explicit ParseGffFasta();
  ~ParseGffFasta();

  // Functionality passed to the implmentation.
  [[nodiscard]] std::shared_ptr<GenomeDatabase> readFastaFile(const std::string& organism, const std::string& fasta_file_name);

  [[nodiscard]] bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);

  [[nodiscard]] bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);

  [[nodiscard]] std::shared_ptr<GenomeDatabase> readFastaGffFile( const std::string& organism,
                                                                  const std::string& fasta_file_name,
                                                                  const std::string& gff_file_name);

  void readTssGffFile(const std::string& tss_gff_file_name, GenomeDatabase& genome_db_ptr);

private:

  class GffFastaImpl;       // Forward declaration of the Gff Fasta parser implementation class
  std::unique_ptr<GffFastaImpl> gff_fasta_impl_ptr_;    // Gff Fasta parser PIMPL

};


}   // end namespace

#endif //KGL_GFF_FASTA_H
