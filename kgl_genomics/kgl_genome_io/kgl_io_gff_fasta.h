///
// Created by kellerberrin on 3/10/17.
//

#ifndef KGL_GFF_FASTA_H
#define KGL_GFF_FASTA_H


#include "kgl_genome_genome.h"
#include "kgl_io_gff3.h"
#include "kgl_io_fasta.h"

#include <memory>
#include <string>


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Reads a fasta and gff file and composes the genome database.
// Static object to provide data hiding and namespace.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseGffFasta {

public:

  ParseGffFasta() = delete;
  ~ParseGffFasta() = delete;

  [[nodiscard]] static std::shared_ptr<GenomeReference> readFastaGffFile( const std::string& organism,
                                                                          const std::string& fasta_file_name,
                                                                          const std::string& gff_file_name);

};


}   // end namespace


#endif
