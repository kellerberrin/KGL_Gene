//
// Created by kellerberrin on 23/08/22.
//

#ifndef KGL_KGL_GFF_FASTA_SEQAN_H
#define KGL_KGL_GFF_FASTA_SEQAN_H

#include "kgl_gff_fasta.h"

namespace kellerberrin::genome {   //  organization::project level namespace



class ParseGffFastaSeqan {

public:

  explicit ParseGffFastaSeqan();
  ~ParseGffFastaSeqan();

  // Functionality passed to the implmentation.
  [[nodiscard]] std::shared_ptr<GenomeReference> readFastaFile(const std::string& organism, const std::string& fasta_file_name);

  [[nodiscard]] bool writeFastaFile(const std::string& fasta_file_name, const std::vector<WriteFastaSequence>& fasta_sequences);

  [[nodiscard]] bool readFastaFile(const std::string& fasta_file_name, std::vector<ReadFastaSequence>& fasta_sequences);

  [[nodiscard]] std::shared_ptr<GenomeReference> readFastaGffFile(const std::string& organism,
                                                                  const std::string& fasta_file_name,
                                                                  const std::string& gff_file_name);

private:

  class GffFastaSeqan;       // Forward declaration of the Gff Fasta parser implementation class

  using GffFastaImpl = GffFastaSeqan;
  std::unique_ptr<GffFastaImpl> gff_fasta_impl_ptr_;    // Gff Fasta parser PIMPL

};



} // end of namespace.

#endif //KGL_KGL_GFF_FASTA_SEQAN_H
