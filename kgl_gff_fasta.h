///
// Created by kellerberrin on 3/10/17.
//

#ifndef KGL_GFF_FASTA_H
#define KGL_GFF_FASTA_H


#include <memory>
#include <string>
#include <map>
#include "kgl_genome_db.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

// Creates an instance of a Genome database object.
// Parses the input gff(3) file and annotates it with a fasta sequence
// This class is a facade; the file reading details are handled by a 3rd party library.
class ParseGffFasta {

public:

  explicit ParseGffFasta(Logger& logger);
  ~ParseGffFasta();

  // Functionality passed to the implmentation.
  std::shared_ptr<GenomeDatabase> readFastaFile(const std::string& fasta_file_name);

  std::shared_ptr<GenomeDatabase> readFastaGffFile(const std::string& fasta_file_name,
                                                   const std::string& gff_file_name);


private:

  Logger& log; // Should be declared first.

  class GffFastaImpl;       // Forward declaration of the Gff Fasta parser implementation class
  std::unique_ptr<GffFastaImpl> gff_fasta_impl_ptr_;    // Gff Fasta parser PIMPL

};


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_GFF_FASTA_H
