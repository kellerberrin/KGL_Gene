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
#include "kgl_logging.h"


namespace kellerberrin::genome {   //  organization::project level namespace

// Creates an instance of a Genome database object.
// Parses the input gff(3) file and annotates it with a fasta sequence
// This class is a facade; the file reading details are handled by a 3rd party library (Seqan).

// Passed in a vector to write a fasta sequence.
// WriteFastaSequence.first is the sequence name.
// WriteFastaSequence.second is the virtual sequence (DNA or Amino).
using WriteFastaSequence = std::pair<std::string, std::shared_ptr<VirtualSequence>>;
// Passed in a vector from fasta sequence read.
// ReadFastaSequence.first is the sequence name.
// ReadFastaSequence.second is a std::string of the fasta sequence.
// Note that AlphabetSequence and AlphabetSequence_t are very different data types.
using ReadFastaSequence = std::pair<std::string, std::shared_ptr<AlphabetSequence_t>>;


class ParseGffFasta {

public:

  explicit ParseGffFasta();
  ~ParseGffFasta();

  // Functionality passed to the implmentation.
  std::shared_ptr<GenomeDatabase> readFastaFile(const std::string& organism, const std::string& fasta_file_name);

  bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);

  bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);

  std::shared_ptr<GenomeDatabase> readFastaGffFile(const std::string& organism,
                                                   const std::string& fasta_file_name,
                                                   const std::string& gff_file_name);

  void readTssGffFile(const std::string& tss_gff_file_name, GenomeDatabase& genome_db_ptr);

private:

  class GffFastaImpl;       // Forward declaration of the Gff Fasta parser implementation class
  std::unique_ptr<GffFastaImpl> gff_fasta_impl_ptr_;    // Gff Fasta parser PIMPL

};


}   // end namespace

#endif //KGL_GFF_FASTA_H
